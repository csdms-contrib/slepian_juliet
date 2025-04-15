function varargout=covgammiosl(th,params,method,ifinv);
% [covg,kgrads]=COVGAMMIOSL(th,params,method,ifinv);
%
% Calculates the covariance of the SCORE, the first derivative of the
% two-dimensional debiased spatial Whittle Matern log-likelihood, including
% blurring and wavenumber correlation required for calculation of estimator
% covariance (see eq. A54, and compare to eq. 139 of Simons & Olhede 2013, and
% see eq. 37 of Guillaumin et al. 2022). A central ingredient is the covariance
% of the (modified) PERIODOGRAM, and calculations are implemented as:
%
% (1) repeated SAMPLING of periodograms from newly simulated data (parallel);
% (2) exact, direct evaluation using DFTMTX method (fast but memory-intensive);
% (3) exact, direct autocovariance DIAGONALS summation (parallel and GPU array).
%
% INPUT:
%
% th        The Matern parameters, [s2 nu rh]
% params    The parameter structure of the grid
% method    1   SAMPLING - approximation simulating realizations for TH from which
%               we form many empirical periodograms; method can be run in parallel;
%           2   DFTMTX - exact calculation of the covariance of the score by
%               applying Isserlis' rule to express the covariance of the 
%               empirical periodogram in terms of the DFT of the autocovariance,
%               requiring the formation of the full cross-covariance matrix and
%               its DFT in memory (NyNx by NyNx);
%           3   DIAGONALS - the same calculation as method 2 but implemented to
%               avoid forming matrices in memory by summing over the diagonals
%               through intensive looping; method can be run in parallel and 
%               the calculation of the diagonal offset built into 
%               BLUROSY/SPATMAT ('efd' for tsto~=[0 0]) are converted
%               into GPU arrays when appropriate and useful.
% ifinv         Indicates which Matern parameters are inverted for and are required in
%               the covariance calculation as an array of three logical flags.
%
% OUTPUT:
%
% covg      The covariance of the score (i.e., the gradient of the blurred
%           debiased Whittle likelihood). Matrix size depends on ifinv, with the 
%           order returned in the general case of ifinv=[1 1 1] being 
%           [s2s2 s2nu s2rh; nus2 nunu nurh; rhs2 rhnu rhrh].
% kgrads    Specific only to method 1: the score calculated for each of the
%           sampled empirical periodograms; used for demos/validation only.
%
% SEE ALSO:
%
% COVTHOSL, MATERNOSY, BLUROSY
%
% EXAMPLES:
% 
% th=[1 0.5 2]; params=[]; params.NyNx=[8 8]; params.dydx=[1 1];
% params.blurs=-1; params.taper=1; ifinv=[1 1 1];
% tic;  covg1=covgammiosl(th,params,1,ifinv); toc
% tic;  covg2       =covgammiosl(th,params,2,ifinv); toc
% tic;  covg3       =covgammiosl(th,params,3,ifinv); toc
% difer(covg2-covg3)
% difer(covg1-covg3,1)
%
% % For studying random deletions and other non-unit tapers
%
% th=[1 0.5 2];params=[];params.NyNx=[16 16];params.dydx=[1 1];
% params.blurs=-1;ifinv=[1 1 1];tpr=ones(params.NyNx);tpr(2,3)=0;
% tpr=ones(params.NyNx);tpr(randi(prod(params.NyNx),[5 1]))=deal(0);
% params.taper=tpr;
% covg1 = covgammiosl(th,params,1,ifinv);
% covg2 = covgammiosl(th,params,2,ifinv);
% covg3 = covgammiosl(th,params,3,ifinv);
% difer(covg2-covg3)
% difer(covg1-covg3,1)
%
% % Let's copy this into BLUROSY as a demo:
% covgammiosl('demo0')
% % Validation of the covariance methods for the 3-parameter Matern case as
% % error ellipses and cross-correlation heatmaps:
% covgammiosl('demo1') 
% % Validation of the covariance methods for the 2-parameter Matern case as
% % error ellipses and cross-correlation heatmaps:
% covgammiosl('demo2') 
% % Demonstration of convergence of SAMPLE method for increasing number of
% % periodograms:
% covgammiosl('demo3')
% % Demonstration of constructing the DFTMTX terms for very small grids:
% covgammiosl('demo4')
% % [IN PROGRESS] Demonstration of indexing for the DIAGONALS method:
% covgammiosl('demo5')
% % Study of the dependence of parameter covariance on grid-size (NyNx):
% covgammiosl('demo6')
% % Study of the dependence of parameter covariance on density (dydx):
% covgammiosl('demo7')
% % Study of the dependence of parameter covariance on missingness (taper),
% % either as random deletions or a growing masked region:
% covgammiosl('demo8')
%
% Last modified by fjsimons-at-princeton.edu, 04/08/2025
% Last modified by olwalbert-at-princeton.edu, 04/08/2025

if ~isstr(th)
  % Set default for ifinv to assume that the covariance of the score for all
  % three Matern parameters is desired
  if ~exist('ifinv','var'); ifinv=[1 1 1]; end

  % Set the calculation method default to SAMPLING so we don't accidentally
  % get into trouble with large fields (costly in memory for method 2, costly
  % in calculation time for method 3)
  if ~exist('method','var');method=1;end

  % The number of Matern parameters
  np=3;
  % To save on calculation time, disable extra verification checks in BLUROSY
  xver=0;

  % Capture the mask and work it into calculations
  if isfield(params,'taper') & numel(params.taper)>1
    % Our field is masked
    Tx=params.taper;
  else
    % Our field is not masked
    if method==2
      % We will calculate the autocorrelation sequence here and will need an
      % array
      Tx=ones(params.NyNx);
    else
      % Otherwise, it is more advantageous to have as a scalar
      Tx=1;
    end
  end
  % The size of the field
  NyNx=params.NyNx;

  % Calculating the score requires the expected blurred periodogram, which we
  % form from BLUROSY through MATERNOSP
  Sbar=maternosp(th,params);

  % We also need blurred partial derivatives of the expected periodogram, which
  % we form from BLUROSY through MATERNOSY; initialize as [prod(NyNx) np] for
  % storage of the partials as vectors for each of the np parameters
  dSbardth=zeros([size(Sbar) np]);
  % One could use MAOSL but that is just using BLUROSY anyway
  for thnd=1:np
    if ifinv(thnd)
        % Only calculate partials for the parameters requested by ifinv
        dSbardthnd=blurosy(th,params,xver,'efs',[],thnd);
        % Store as vectors
        dSbardth(:,:,thnd)=dSbardthnd;
    end
  end

  % Calculate the covariance of the score for the selected method
  switch method
    case 1
      % Calculate the covariance of the score by SAMPLING many realizations of
      % the empirical periodogram

      % Set the number of realizations (1e3-2e3 is sufficient for terms to
      % asymptote out; if not runinpar, 1e3 takes 8.4s for a 128 by 128 grid, 
      % if runinpar with 23 workers, 1e3 takes 1.1s, with more scaling linearly)
      numreals=1000; 
      % Pre-allocate space for each score
      kgrads=zeros(numreals,np);

      % Simulate repetitively for the same parameters using circulant embedding,
      % returning the empirical periodogram of each realization, Semp
      runinpar=0;
      % Simulate via circulant embedding
      simparams=params; simparams.blurs=Inf;
      if runinpar==1
        % Create the parallel pool with NumWorkers; running into problems using
        % all physical cores, despite twice as many logical cores detected
        NumWorkers=feature('numcores')-1;
        if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end
        parfor ind=1:numreals
          % Generate an empirical (modified, tapered, windowed) periodogram
          % through SGP simulation; SIMULOSL already takes the taper in to
          % account and normalizes Semp appropriately
          [~,~,~,~,~,Semp]=simulosl(th,simparams,xver);
          Semps(:,ind)=Semp;
        end
      else
          % Same as above, but without PARFOR
        Semps=zeros(numel(Sbar),numreals);
        for ind=1:numreals
          [Hx,~,~,k,Hk,Semp]=simulosl(th,simparams,xver);
          Semps(:,ind)=Semp;
          % Technically could use GAMMIOSL, if we fed it k and Hk
          ggrads(:,ind)=gammiosl(k,th,params,Hk,xver);
        end
        ggrads=ggrads';
      end
      
      % Calculate the elements of the score for each parameter vectorized in
      % dSbardth [numreals by np] % remember dSbardth./Sbar.^2 <- MAOSL/Sbar
      % (A.53) with Sbar factored out; (A.54) would have gone via cov(semps)
      % kgrads=deal(squeeze(sum((1-Semps/Sbar).*mths)));
      % kgrads=deal(squeeze(sum((1-Semps/Sbar).*dSbardth./Sbar)));
      kpp=knums(params);
      % Take out zero wavenumber because we only ever  simulate
      % ZERO-MEAN PROCEESSES and only ever analyze DEMEANED data
      kgrads=deal(squeeze(sum((Sbar(~~kpp)-Semps(~~kpp,:)).*...
                              dSbardth(~~kpp,:,:)./Sbar(~~kpp).^2))); 
      lk=length(kpp(~~kpp));

      % The fair comparison would be between, essentially identical
      % diferm(ggrads,kgrads/lk)
      % [cov(ggrads) ; cov(kgrads/lk]
      % But in GAMMIOSL the zero wavenumber was excluded but in here not
      
      % Normalization factors in the case of masked and unit tapers
      if isfield(params,'taper') & numel(params.taper)==prod(params.NyNx)
         % Adjust for the size of an explicit taper as in SIMULOSL and MLEOSL
         % OLW/FJS Need to think about implications of removing a zero wavenumber
         normfact=(prod(params.NyNx)./sum(Tx(:).^2)).^2.*prod(params.NyNx).^2;
      else
         normfact=lk^2;
      end
      % Straight calculation of the covariance of the score
      covg=cov(kgrads)./normfact;
    case 2
      % Exact calculation of the gradient of the empirical periodograms using
      % Isserlis' Rule, which provides an equivalence with the linear combination
      % of transposed and conjugated Fourier Transforms of the autocovariance
      % (see Walden+1999 eqs 7-8):
      %   In 1-D, cov{|Hk|^2,|Hk'|^2} = U'*Cmn*U + U'*Cmn*U.'
      % We calculate the Matern autocovariance calculated a NyNx^2 by NyNx^2
      % distance grid that incorporates all possible lag pairs, which in the 
      % frequency domain represents all finite, discrete wavenumber interactions
      % for the data observation grid. 
      % Note that all matrices are completely assembled in memory using this 
      % method; the largest grid size I can calculate using R2024a on Aguilonius 
      % is about 128x128, anything larger is killed. Method 3 works equivalently 
      % and is effective for larger grid sizes, but is a slower calculation.
      % In the future, we might think about whether we can avoid forming all in
      % memory and taking advantage of the symmetry of the dftmtx and the
      % autocovariance matrix: Cmn is symmetric and positive-definite, which
      % allows for Cholesky decomposition
      % chCmn=chol(Cmn);

      % The distance grid, dxxp, of all pairwise L2 distances between two common 
      % gridded periodograms assembled as a block Toeplitz matrix, taking into 
      % account the spacing in x and y
      dydx=params.dydx;NyNx=params.NyNx;
      ys=0:NyNx(1)-1;xs=0:NyNx(2)-1;
      [grdx,grdy]=meshgrid(xs(:),ys(:));
      grdx=grdx(:).*dydx(2);
      grdy=grdy(:).*dydx(1);
      lagx=grdx-grdx';
      lagy=grdy-grdy';
      dxxp=sqrt(lagx.^2+lagy.^2);
      % XXPDIST provides an equivalent distance matrix:
      % dxxp=xxpdist([grdx grdy]);

      % Calculate the autocovariance at the lag-distances for the Matern model
      Cmn=maternosy(dxxp,th);
      
      % If we provided a non-unit taper, we must incorporate it in our analysis 
      % by taking its product with our spatial autocovariance matrix
      if numel(params.taper)>1                                                  
        % This is what I thought the solution should have been:
        % Reconstruct the mask, Tx, so that it corresponds with the entries of the
        % autocovariance matrix, Cmn
        % for ind=1:NyNx(1)                                                       
        %   begi=1+(ind-1)*NyNx(1);                                               
        %   endi=NyNx(1)+(ind-1)*NyNx(1);                                         
        %   tb=toeplitz(Tx(begi:endi));                                           
        %   nxnxcircblock{ind}=tb;
        % end                                                                     
        % tt = cell2mat(nxnxcircblock(toeplitz(1:NyNx(1))));                      
        % Cmn=Cmn.*tt;                                                          
        % This doesn't always work:
        % Cmn=matslice(matslice(ones(size(Cmn)),Tx(:)),Tx(:),-1).*Cmn;          
        % This is what matches method 3:
        Cmn=Cmn.*Tx(:).*Tx(:)'; 
      end

      % In 1-D, for a vector m x 1
      %   fft(Cm)  =(Cm.'*Um).' + O(eps)
      % In 2-D, for an m x n matrix
      %   fft2(Cmn)=Um*Cmn*Un.' + O(eps)
      %            =tensorprod(Un,tensorprod(Um,Cmn.',2,2).',2,1).'+O(eps)
      %            =tensorprod(Un,tensorprod(Um,Cmn  ,2,1)  ,2,2).'+O(eps)
      % The first term of the covariance of the periodogram via Isserlis' rule 
      % (eq. 9 of VoWE; also see Walden+1994)
      % .' is TRANSPOSE whereas ' is CTRANSPOSE
      % In 1-D, EJJht=Um'*Cm*Um./n;
      % In 2-D, EJJht=(Um*(Um*Cmn*Un.')'*Un.')'/prod(NyNx);
      % We have Ny*Nx sets of autocovariance blocks that we will have to reorder
      % correctly after each operation

      % The distance matrix is comprised of Nx [Ny by Ny] blocks in the rows and
      % columns. Reshape Cmn so that there are NyNx ``column vectors'' that are 
      % Nx sets of the [Ny by Ny] blocks; dimensions should be [Nx by Ny by NyNx]
      Cm_n_mn=reshape(Cmn,[NyNx(1) NyNx(2) prod(NyNx)]);
      % The DFT matrix (mxm)
      Um=dftmtx(size(Cm_n_mn,1));
      % The DFT matrix (nxn)
      Un=dftmtx(size(Cm_n_mn,2));
      % Form the shared inner term of both EJJht and EJJt
      % This is like taking fft2(Cm_n_mn(:,:,i),[],2) for the ith slice
      EJJht_in=reshape(permute(...                                              
                tensorprod(Um,tensorprod(Un,Cm_n_mn,1,2),1,2),...               
               [3 1 2]),NyNx(1),NyNx(2),prod(NyNx)); 

      % Form the outer term of EJJht
      % SECOND term in right hand side of A57 "Hermitian  transpose"
      EJJht_out=conj(tensorprod(Um,tensorprod(Un,conj(EJJht_in),1,2),1,2));
      EJJht_out=reshape(permute(EJJht_out,[3 1 2]),[prod(NyNx) prod(NyNx)]);
      diferm(EJJht_out-EJJht_out.')

      % The second term of the covariance of the periodogram via Isserlis' rule 
      % (eq. 9 of VoWE; also see Walden+1994)
      % In 1-D, EJJt=Um'*Cm*Um.'./NyNx;
      % In 2-D, EJJt=(Um*(Um*Cmn*Un.')*Un)'/prod(NyNx)
      % FIRST term in right hand side of A57 "transpose"
      EJJt_out=tensorprod(Um,tensorprod(Un,EJJht_in,1,2),1,2);
      EJJt_out=reshape(permute(EJJt_out,[3 1 2]),[prod(NyNx) prod(NyNx)]);
      diferm(EJJht_out-EJJht_out') 
      
      % Sum and normalize the quantities above; eq. 8 of Walden+1994
      Hk2cov=(abs(EJJht_out).^2+abs(EJJt_out).^2)./sum(Tx,"all").^2;

      % Un-do the normalization factor applied in BLUROSY for the blurred
      % expected periodogram and its partial derivatives; shift both such
      % that the zero-wavenumber is in the (1,1) position 
      nf=(2*pi)^2;
      Sbar=ifftshift(v2s(Sbar,params))*nf;
      % [NyNx by 1]
      Sbar=Sbar(:);
      for ind=1:np
        tmp=ifftshift(v2s(dSbardth(:,:,ind),params));
        dSbardth(:,:,ind)=tmp(:)*nf; 
      end
      % [NyNx by np]
      dSbardth=squeeze(dSbardth);

      % Partial derivative divided by square of periodogram, mthth/Sbar
      fax=dSbardth./Sbar.^2;

      % Optional: normalize by another Sbar.^2; Arthur special ?
      normalize=0;
      if normalize
        Sbar_n=v2s(Sbar,params);
        seq_ep=Sbar_n.*Sbar_n;
      else
        seq_ep=1;
      end 

      % Normalization factors for both non-unit and unit tapers
      if isfield(params,'taper')&numel(params.taper)==prod(NyNx)
        normfact=(prod(params.NyNx)./sum(Tx(:).^2)).^2.*prod(params.NyNx).^2;
      else
        normfact=prod(NyNx).^2;
      end

      % Final assembly (A54)
      covg=fax'*Hk2cov*fax./seq_ep./normfact;
      % I know that I need to multiply by prod(dydx).^2, but I am not sure where
      % would be best. Doing so here works for now.
      covg=covg.*prod(dydx).^2;
      % Uncomment below to see  the number of bytes each variable occupies in
      % the workspace
      % wkspc=whos;sum([wkspc.bytes])
    case 3
      % Calculates the covariance of the score using the per-diagonal FFT method
      % Calc times after improving efficiency: 48x48 0.7s; 128x128 33s; 256x256 970s... 

      % Currently developed following the organization and normalization
      % strategy that Arthur uses, meaning that we will multiply our
      % periodograms and their derivatives by (2*pi)^2, shift them so that 0 lag
      % and wavenumber is in the C11 position rather than centered, and we will 
      % call the new BLUROSY method 'efd' to allow for offsets in wavenumber
      nf=(2*pi)^2;
      Sbar=ifftshift(v2s(Sbar,params))*nf;
      Sbar=Sbar(:);
      for ind=1:np
        tmp=ifftshift(v2s(dSbardth(:,:,ind),params));
        dSbardth(:,:,ind)=tmp(:)*nf; 
      end
      dSbardth=squeeze(dSbardth);
      % Partial derivative divided by square of periodogram, mthth/Sbar
      fax=dSbardth./Sbar.^2;

      % To offset the periodograms and their derivatives, we calculate the index
      % offset based on a NyNx grid
      Ny=NyNx(1);
      Nx=NyNx(2);

      % Ask Arthur when we should (not) normalize; default is false
      normalize=0;
      if normalize
        Sbar_n=v2s(Sbar,params);
      else
        Sbar_n=ones(NyNx);
        seq_ep =1;
        seq_ep2=1;
      end
      Sbar_n=Sbar_n(:);

      % Create a 0...Nx-1 by 0...Ny-1 meshgrid for indexing the offset of the
      % FFT of the spatial autocovariance 
      [x,y]=meshgrid(0:Nx-1,0:Ny-1);
      % Initialize array of summed diagonals
      s_1=0;
      s_2=0;

      % Run in parallel?
      runinpar=1;
      if runinpar
        % Create the parallel pool with NumWorkers; should be 23 on Aguilonius
        NumWorkers=feature('numcores')-1;
        if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end
        parfor my=(-Ny+1:Ny-1)
          for mx=(-Nx+1:Nx-1)                                                    
            % Create two logical vectors for indexing the offset in diagonals
            % of the FFT of the autocovariance
            a=y<Ny-my&x<Nx-mx;                                                
            b=y>=-my&x>=-mx;                                                  
            ep_ind1=a&b;                                                      
            a=y>=my&x>=mx;                                                    
            b=y<Ny+my&x<Nx+mx;                                                
            ep_ind2=a&b;                                                      

            % Create a vector of logicals to index for offset in the 
            % anti-diagonals of the FFT of the autocovariance
            mx2=mx+Nx-1;
            my2=my+Ny-1;
            a=x<=mx2&y<=my2;
            b=x>=mx2-Nx+1&y>=my2-Ny+1;
            ep_ind3=a&b;
            ep_ind3=ep_ind3(:);

            % We acquire the offset FFT of the autocovariance via BLUROSY, which
            % we now provide method 'efd' and indices to calculate the offsets
            % in frequency as the tsto input argument
            fftcov_diag=blurosy(th,params,xver,'efd',[my mx]);                        
            fftcov_diag=v2s(fftcov_diag,params)*nf;                                           
            cdd=fftcov_diag(max([1,my+1]):min([Ny+my,Ny]),...                         
                            max([1,mx+1]):min([Nx+mx,Nx]));                           
            cdd_l2=abs(cdd).^2;                                                  
            cdd_l2=cdd_l2(:);                                                       
                                                                                 
            % Use the same framework to select the antidiagonals from BLUROSY
            fftcov_anti=blurosy(th,params,xver,'efd',[my2 mx2]);
            fftcov_anti=v2s(fftcov_anti,params)*nf;
            cad=fftcov_anti(max([1,my2-Ny+2]):min([my2+1,Ny]),...
                            max([1,mx2-Nx+2]):min([mx2+1,Nx]));
            cad_l2=abs(cad).^2;
            cad_l2=cad_l2(:);

            % We need to be a bit more verbose about this in PARFOR
            if normalize
              seq_ep=Sbar_n(ep_ind1).*Sbar_n(ep_ind2);
              Sbar_n3=Sbar_n(ep_ind3);
              seq_ep2=Sbar_n3.*flip(Sbar_n3,1);
            else
              seq_ep =1;
              seq_ep2=1;
            end

            % The partials indexed for the current `diagonal'
            f1=fax(ep_ind1,:); f2=fax(ep_ind2,:);
            % Product of all partials in row-order:
            % [s2s2 s2nu s2rh nus2 nunu nurh rhs2 rhnu rhrh]
            f_seq=[f1(:,1).*f2 f1(:,2).*f2 f1(:,3).*f2];
            s_1=s_1+sum(cdd_l2.*f_seq./seq_ep,1);

            % Arthur tends to comment this out and double s_1 instead.
            % Are s_1 and s_2 always the same? Not quite, but they probably sum
            % to the same thing in the end.
            f1=fax(ep_ind3,:);
            f_seq2=[f1(:,1).*flipud(f1) f1(:,2).*flipud(f1) f1(:,3).*flipud(f1)];
            s_2=s_2+sum(cad_l2.*f_seq2./seq_ep2,1);
          end
        end
      else
          % Not a parallel run
          % Useful for debugging; calculation time depends on grid-size
          % Display progress?
          dispprog=0;
          rnd=0;
          for my=(-Ny+1:Ny-1)
              if dispprog
                  disp(sprintf('%i%% complete',...
                               floor(((my+Ny-1)/(2*Ny)*100))))
              end
              for mx=(-Nx+1:Nx-1)                                                    
                  rnd=rnd+1;
                  % Create two logical vectors for indexing the offset in diagonals
                  % of the FFT of the autocovariance
                  a=y<Ny-my&x<Nx-mx;                                                
                  b=y>=-my&x>=-mx;                                                  
                  ep_ind1=a&b; ep_ind1=ep_ind1(:);
                  a=y>=my&x>=mx;                                                    
                  b=y<Ny+my&x<Nx+mx;                                                
                  ep_ind2=a&b; ep_ind2=ep_ind2(:);

                  % Create a vector of logicals to index for offset in the 
                  % anti-diagonals of the FFT of the autocovariance
                  mx2=mx+Nx-1;
                  my2=my+Ny-1;
                  a=x<=mx2&y<=my2;
                  b=x>=mx2-Nx+1&y>=my2-Ny+1;
                  ep_ind3=a&b; ep_ind3=ep_ind3(:);

                  % We acquire the offset FFT of the autocovariance via BLUROSY, which
                  % we now provide method 'efd' and indices to calculate the offsets
                  % in frequency as the tsto input argument
                  fftcov_diag=blurosy(th,params,xver,'efd',[my mx]);                        
                  fftcov_diag=v2s(fftcov_diag,params)*nf;                                           
                  cdd{rnd}=fftcov_diag(max([1,my+1]):min([Ny+my,Ny]),...                         
                                  max([1,mx+1]):min([Nx+mx,Nx]))
                  cdd_l2=abs(cdd{rnd}).^2;                                                  
                  cdd_l2=cdd_l2(:);                                                       
                  
                  % Use the same framework to select the antidiagonals from BLUROSY
                  fftcov_anti=blurosy(th,params,xver,'efd',[my2 mx2]);
                  fftcov_anti=v2s(fftcov_anti,params)*nf;
                  cad{rnd}=fftcov_anti(max([1,my2-Ny+2]):min([my2+1,Ny]),...
                                  max([1,mx2-Nx+2]):min([mx2+1,Nx]))
                  cad_l2=abs(cad{rnd}).^2;
                  cad_l2=cad_l2(:);

                  if normalize
                      seq_ep =Sbar_n(ep_ind1).*Sbar_n(ep_ind2);
                      ep_n3=ep_n(ep_ind3);
                      seq_ep2=ep_n3.*flip(ep_n3,1);
                  end

                  % The partials indexed for the current `diagonal'
                  f1=fax(ep_ind1,:); f2=fax(ep_ind2,:);
                  % Product of all partials in row-order:
                  % [s2s2 s2nu s2rh nus2 nunu nurh rhs2 rhnu rhrh]
                  f_seq=[f1(:,1).*f2 f1(:,2).*f2 f1(:,3).*f2];
                  s_1=s_1+sum(cdd_l2.*f_seq./seq_ep,1);

                  % Arthur tends to comment this out and double s_1 instead. Are s_1
                  % and s_2 always the same?
                  f1=fax(ep_ind3,:);
                  f_seq2=[f1(:,1).*flipud(f1) f1(:,2).*flipud(f1) f1(:,3).*flipud(f1)];
                  s_2=s_2+sum(cad_l2.*f_seq2./seq_ep2,1);

                  % We want to study how different s_1 and s_2 are; we see that stepping
                  % through the offset periodogram yields very different values, which
                  % is expected due to different indices (ep_ind1 & ep_ind2 vs
                  % ep_ind3), but we are really interested in whether s_1 and s_2 are
                  % similar following the complete summation.
                  % disp(sprintf('s_1 - s_2: %g %g %g %g %g %g %g %g %g',s_1-s_2))
              end
          end
      end
      % for ind=1:floor(size(cdd,2)/2); difer(cdd{ind}-conj(cdd{end-ind+1}));
      % end

     % Following several tests, we see that there is a very small difference 
     % between s_1 and s_2 following summation (O(-10) - O(-16)), with
     % dependence on th. Probably nothing to keep worrying about, and no gains expected
     % whether you deal with it or not
     % disp(sprintf('s_1 - s_2: %g %g %g %g %g %g %g %g %g',s_1-s_2))
     s=reshape(s_1+s_2,np,np);

     % Normalization factors in the case of masked and unit tapers
     if isfield(params,'taper') & numel(params.taper)==prod(NyNx)
        normfact=(prod(params.NyNx)./sum(Tx(:).^2)).^2.*prod(params.NyNx).^2;
     else
         normfact=prod(params.NyNx).^2;
     end
     covg=s./normfact;
     % wkspc=whos; sum([wkspc.bytes])
  end % end of method cases
  % Return covg only for the parameters that were involved in the inversion
  covg=matslice(covg,ifinv);

  % If we used method 2 or 3, grads was never defined
  defval('grads',NaN)
  % Optional output
  varns={covg,grads};
  varargout=varns(1:nargout);
elseif strcmp(th,'demo0')
    % SHOULD THIS BE A BLUROSY DEMO?
    % Creates a plot of an empirical periodogram, blurred expected periodogram, 
    % the spatial covariance, and their partial derivatives
    th=[1 0.5 2];params=[];params.NyNx=[32 32];params.dydx=[1 1];
    params.blurs=Inf;
    [~,~,~,~,Hk]=simulosl(th,params,0);
    params.blurs=-1;
    [Sbar,~,~,acv]=blurosy(th,params);
    [dSbards2,~,~,d_acv_s2]=blurosy(th,params,[],[],[],1);
    [dSbardnu,~,~,d_acv_nu]=blurosy(th,params,[],[],[],2);
    [dSbardrh,~,~,d_acv_rh]=blurosy(th,params,[],[],[],3);
    set(groot, 'DefaultTextInterpreter', 'latex')
    clf;
    [ah,ha,H]=krijetem(subnum(3,3));
    axes(ha(1));[~,cax]=imagefnan([],[],v2s(abs(Hk).^2),[],[-1 1]);title('$|Hk^2|$')
    axes(ha(2));imagefnan([],[],v2s(Sbar),[],cax);title('$\bar{S}(k)$')
    axes(ha(3));imagefnan([],[],v2s(acv),[],cax);title('$C(y)$')
    axes(ha(4));imagefnan([],[],v2s(dSbards2),[],cax);title('$d\bar{S}(k)/d\sigma^2$');
    axes(ha(5));imagefnan([],[],v2s(dSbardnu),[],cax);title('$d\bar{S}(k)/d\nu$');
    axes(ha(6));imagefnan([],[],v2s(dSbardrh),[],cax);title('$d\bar{S}(k)/d\rho$');
    axes(ha(7));imagefnan([],[],v2s(d_acv_s2),[],cax);title('$dC(y)/d\sigma^2$');
    axes(ha(8));imagefnan([],[],v2s(d_acv_nu),[],cax);title('$dC(y)/d\nu$');
    axes(ha(9));imagefnan([],[],v2s(d_acv_rh),[],cax);title('$dC(y)/d\rho$');
    for ind=1:9
      axes(ah(ind))
      longticks
    end
    cb=colorbar;
    cb.Limits=cax; cb.Location="southoutside"; cb.Position=[0.35 0.035 0.3 0.01];
    movev(ha,-0.025) 
    ti=sgtitle(sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s',...
       '\sigma^2=',th(1),' \mathrm{[unit]}^2, \nu=',th(2),', \rho=',th(3),...
       ' m; ',params.NyNx(1).*params.dydx(1),' m x ',...
       params.NyNx(2).*params.dydx(2),'m'),... 
       'Interpreter','latex');
elseif strcmp(th,'demo1')
   % Evaluating COVGAMMIOSL for inversion of all three Matern parameters as
   % cross-plots with error ellipses for the four methods of calculating
   % parameter covariance, and as correlation plots 
   ifinv=[1 1 1]; 

   try
     % Load a previous simulation suite
     datum='18-Feb-2025';
     % mleosl('demo2',datum)
     [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(datum,100);
     if isinf(p.blurs); p.blurs=-1; end
     if ~isfield(p,'taper'); p.taper=1; end
   catch
     % Or calculate a new one
     p=[];p.NyNx=[56 56];p.dydx=[1 1];p.blurs=-1;p.taper=1;
     th0=[3 1.5 5];
     numreals=300;
     mleosl('demo1',numreals,[3 1.5 5],p)
     datum=date;
     thhats=load(sprintf('mleosl_thhat_%s',datum));
   end

   % Calculate parameter covariance for each method, storing the cross-parameter
   % covariance as unique variables
   covg1=covgammiosl(th0,p,1,ifinv);
   cov1=covthosl(th0,p,covg1,ifinv);
   cov1s2nu=matslice(cov1,[1 1 0]);
   cov1s2rh=matslice(cov1,[1 0 1]);
   cov1nurh=matslice(cov1,[0 1 1]);

   covg2=covgammiosl(th0,p,2,ifinv);
   cov2=covthosl(th0,p,covg2,ifinv);
   cov2s2nu=matslice(cov2,[1 1 0]);
   cov2s2rh=matslice(cov2,[1 0 1]);
   cov2nurh=matslice(cov2,[0 1 1]);

   covg3=covgammiosl(th0,p,3,ifinv);
   cov3=covthosl(th0,p,covg3,ifinv);
   cov3s2nu=matslice(cov3,[1 1 0]);
   cov3s2rh=matslice(cov3,[1 0 1]);
   cov3nurh=matslice(cov3,[0 1 1]);

   cove=cov(thhats);
   coves2nu=matslice(cove,[1 1 0]);
   coves2rh=matslice(cove,[1 0 1]);
   covenurh=matslice(cove,[0 1 1]);

   % Means
   mobs1=nanmean(thhats(:,[1 2]));
   mobs2=nanmean(thhats(:,[1 3]));
   mobs3=nanmean(thhats(:,[2 3]));

   % Set-up the 95% confidence level error-ellipses 
   t=linspace(0,2*pi);
   cl=0.95;
   s=chi2inv(cl,2);

   % Calculate the error ellipses whose semi-axes length and orientations depend
   % on the SVD of the covariance matrices in the 2-parameter projected spaces
   [V,D]=eig(cov1s2nu);
   a11=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(cov1s2rh);
   a12=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(cov1nurh);
   a13=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   
   [V,D]=eig(cov2s2nu);
   a21=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(cov2s2rh);
   a22=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(cov2nurh);
   a23=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   
   [V,D]=eig(cov3s2nu);
   a31=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(cov3s2rh);
   a32=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(cov3nurh);
   a33=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];

   [V,D]=eig(coves2nu);
   ae1=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(coves2rh);
   ae2=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   [V,D]=eig(covenurh);
   ae3=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   
   clf;
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];
   % Plot of s2 vs nu
   subplot(131)
   plot(thhats(:,1),thhats(:,2),'ko')
   hold on
   e1(1)=plot(a11(1,:)+mobs1(1),a11(2,:)+mobs1(2),'Color',clrs(1),'LineWidth',1);
   e2(1)=plot(a21(1,:)+mobs1(1),a21(2,:)+mobs1(2),'Color',clrs(2),'LineWidth',1);
   e3(1)=plot(a31(1,:)+mobs1(1),a31(2,:)+mobs1(2),'Color',clrs(3),'LineWidth',1);
   ee(1)=plot(ae1(1,:)+mobs1(1),ae1(2,:)+mobs1(2),'Color','k','LineWidth',1);
   longticks
   xlabel('variance, $\sigma^2$','Interpreter','latex')
   ylabel('smoothness, $\nu$','Interpreter','latex')
   % Plot of s2 vs rho 
   subplot(132)
   plot(thhats(:,1),thhats(:,3),'ko')
   hold on
   e1(2)=plot(a12(1,:)+mobs2(1),a12(2,:)+mobs2(2),'Color',clrs(1),'LineWidth',1);
   e2(2)=plot(a22(1,:)+mobs2(1),a22(2,:)+mobs2(2),'Color',clrs(2),'LineWidth',1);
   e3(2)=plot(a32(1,:)+mobs2(1),a32(2,:)+mobs2(2),'Color',clrs(3),'LineWidth',1);
   ee(2)=plot(ae2(1,:)+mobs2(1),ae2(2,:)+mobs2(2),'Color','k','LineWidth',1);
   longticks
   xlabel('variance, $\sigma^2$','Interpreter','latex')
   ylabel('range, $\rho$','Interpreter','latex')
   % Plot of nu vs rho 
   subplot(133)
   plot(thhats(:,2),thhats(:,3),'ko')
   hold on
   e1(3)=plot(a13(1,:)+mobs3(1),a13(2,:)+mobs3(2),'Color',clrs(1),'LineWidth',1);
   e2(3)=plot(a23(1,:)+mobs3(1),a23(2,:)+mobs3(2),'Color',clrs(2),'LineWidth',1);
   e3(3)=plot(a33(1,:)+mobs3(1),a33(2,:)+mobs3(2),'Color',clrs(3),'LineWidth',1);
   ee(3)=plot(ae3(1,:)+mobs3(1),ae3(2,:)+mobs3(2),'Color','k','LineWidth',1);
   longticks
   xlabel('smoothness, $\nu$','Interpreter','latex')
   ylabel('range, $\rho$','Interpreter','latex')
   [leg1,legic]=legend(sprintf('%i \\theta',size(thhats,1)),'sample','dftmtx',...
            'diagonals','empirical','box','off')
   % Find the 'line' objects
   legicms = findobj(legic,'Type','line');
   for ind=1:size(legicms,1)
     if size(legicms(ind).XData,1)==1
        legicms(ind).XData=0.3;
     else
        legicms(ind).XData(1)=0.2;
        legicms(ind).XData(2)=0.4;
     end
   end
   % Resize the marker in the legend
   %set(legicls,'MarkerSize',20);
   sti=sgtitle(sprintf('$\sigma^2=%0.2f, \nu=%0.2f, \rho=%0.2f$ m; %im x %im',...
           th0,p.NyNx.*p.dydx),'Interpreter','latex');

   % Pause to confirm the figure is nice before saving
   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo1_3pcrossplots_%ix%i_s2%inu%irh%i_%s.eps',...
       p.NyNx.*p.dydx,th0.*[1 10 1],date),'epsc')
   
   % Plot the normalized covariance heatmaps for the four methods
   figure()
   kelicol
   cmrange=[-1 1];
   labs={'s2','nu','rho'};
   labs=matslice(labs,ifinv);
   np=sum(ifinv);
   % The empirical covariance
   ah(1)=subplot(221);
   im(1)=imagesc(matslice(cove,ifinv)./...
         [diag(sqrt(matslice(cove,ifinv)))*...
          diag(sqrt(matslice(cove,ifinv)))'],...
          cmrange);
   title('empirical')
   axis image
   % The sample covariance
   ah(2)=subplot(222);
   im(2)=imagesc(cov1./[diag(sqrt(cov1))*diag(sqrt(cov1))'],cmrange);
   title('sample')
   axis image
   % The diagonals covariance
   ah(3)=subplot(223);
   im(3)=imagesc(cov2./[diag(sqrt(cov2))*diag(sqrt(cov2))'],cmrange);
   title('diagonals')
   axis image
   % The dftmtx covariance
   ah(4)=subplot(224);
   im(4)=imagesc(cov3./[diag(sqrt(cov3))*diag(sqrt(cov3))'],cmrange);
   title('dftmtx')
   axis image
   % Labels, colorbar, positioning, title
   set(ah(:),'xtick',1:np,'XTickLabel',labs)
   set(ah(:),'ytick',1:np,'YTickLabel',labs)
   longticks(ah(:))
   cb=colorbar;
   cb.Label.String='normalized';
   moveh(ah(2),-0.04)
   movev(ah(1),-0.02)
   ah(2).Position(4)=ah(1).Position(4);ah(2).Position(3)=ah(1).Position(3);
   ah(2).Position(2)=ah(1).Position(2);
   ti=sgtitle(sprintf('\\sigma^2=%0.2f, \\nu=%0.2f, \\rho=%0.2f; %im x %im',...
           th0,p.NyNx.*p.dydx))

   % Pause to confirm the figure is nice before saving
   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo1_3pcorrplots_%ix%i_s2%inu%irh%i_%s.eps',...
       p.NyNx.*p.dydx,th0.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo2')
    % Evaluating COVGAMMIOSL for inversion of only two parameters for direct
    % comparison with diagnostics I created for Arthur's jmatrix_sample. 
    % Outputted figure plots kdes of the covariance marginals over histograms of
    % parameter estimates, error ellipses calculated from the parameter 
    % covariances on a cross-plot of parameter estimates, and the asymptotic
    % behavior of the SAMPLE method for increasing number of simulated periodograms
    
    % We need to acquire an MLEOSL suite, either by loading or fresh calculation
    % th1=[1 1/2 5]; 
    % th1=[10 1/2 10]; 
    th1=[5 1/2 7]; 
    params=[]; params.NyNx=[94 97]; params.dydx=[1 1];
    params.blurs=-1; params.taper=1;
    ifinv=[1 0 1];
    datum=date;
    try 
      thhats=load(sprintf('mleosl_thhat_%s',datum));
    catch
      numreals=200;
      mleosl('demo1',numreals,th1,params,[],[],th1,[1 0 1])
      thhats=load(sprintf('mleosl_thhat_%s',datum));
    end
    covemp=cov(thhats);
    mobs=nanmean(thhats);

    % For the covgammiosl sampling method, we want to study the effect of increasing
    % number of samples and will recalculate the variance using F1inv depending
    % on the number of samples
    simparams=params; simparams.blurs=Inf;
    [~,~,~,k]=simulosl(th1,simparams,1);
    F1=fishiosl(k,th1,params);
    F1=matslice(F1,ifinv);
    covF1=inv(F1);
    [~,covg1grads]=covgammiosl(th1,params,1,ifinv);
    numsims=size(covg1grads,1);
    covg1s2s2=zeros(numsims,1); 
    covg1s2rh=zeros(numsims,1);
    covg1rhrh=zeros(numsims,1);
    cov1s2s2=zeros(numsims,1); 
    cov1s2rh=zeros(numsims,1);
    cov1rhrh=zeros(numsims,1);
    for ind=2:numsims
        covg1=real(cov(covg1grads(1:ind,:)))./prod(params.NyNx).^2;
        covg1=matslice(covg1,[1 0 1]);
        cov1=covF1*covg1*covF1;
        covg1s2s2(ind,:)=covg1(1,1);
        covg1s2rh(ind,:)=covg1(1,end);
        covg1rhrh(ind,:)=covg1(end,end);
        cov1s2s2(ind,:)=cov1(1,1);
        cov1s2rh(ind,:)=cov1(1,end);
        cov1rhrh(ind,:)=cov1(end,end);
    end
    % Calculate the exact variance using the diagonal method (2) and the dftmtx 
    % method (3)
    covg2=covgammiosl(th1,params,2,ifinv);
    covg3=covgammiosl(th1,params,3,ifinv);
    cov2=covthosl(th1,params,covg2,ifinv);
    cov3=covthosl(th1,params,covg3,ifinv);

    % Plot
    clf;
    % Color friendly palette
    clrs=["#1B9E77","#D95F02","#7570B3"];

    % Histogram of s2 estimates
    subplot(221)
    numbins=ceil(log2(size(thhats,1)))+1;
    hpl(1)=histogram(thhats(:,1),numbins,'Normalization','pdf');
    hpl(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    vecmod=[0.95 1.05];
    veclen=numsims*10;
    sigvec=linspace(min(thhats(:,1))*vecmod(1),max(thhats(:,1))*vecmod(2),veclen);
    rhovec=linspace(min(thhats(:,3))*vecmod(1),max(thhats(:,3))*vecmod(2),veclen);
    hold on
    vpl(1)=xline(th1(:,1),'Label','\theta_0','LabelOrientation','horizontal');
    vpl(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    nepl(1)=plot(sigvec,normpdf(sigvec,mobs(:,1),sqrt(covemp(1,1))),...
            'k','LineWidth',2,'DisplayName','empirical');
    nspl(1)=plot(sigvec,normpdf(sigvec,th1(:,1),sqrt(cov1s2s2(end))),...
            'Color',clrs{1},'LineWidth',2,'DisplayName','sample');
    ndpl(1)=plot(sigvec,normpdf(sigvec,th1(:,1),sqrt(cov2(1,1))),':',...
            'Color',clrs{2},'LineWidth',2,'DisplayName','diagonals');
    nfpl(1)=plot(sigvec,normpdf(sigvec,th1(:,1),sqrt(cov3(1,1))),'--',...
            'Color',clrs{3},'LineWidth',2,'DisplayName','dftmtx');
    xlabel('$\sigma^2$')
    ylabel('density')
    lg1=legend()
    lg1.AutoUpdate='off';

    % Histogram of rh estimates
    subplot(222)
    hpl(2)=histogram(thhats(:,3),numbins,'Normalization','pdf');
    vpl(2)=xline(th1(:,3),'Label','\theta_0','LabelOrientation','horizontal');
    vecmod=[0.95 1.05];
    veclen=numsims*10;
    hold on
    nepl(2)=plot(rhovec,normpdf(rhovec,mobs(:,3),sqrt(covemp(3,3))),...
                'Color','k','LineWidth',2);
    nspl(2)=plot(rhovec,normpdf(rhovec,th1(:,3),sqrt(cov1rhrh(end))),...
                'Color',clrs{1},'LineWidth',2);
    ndpl(2)=plot(rhovec,normpdf(rhovec,th1(:,3),sqrt(cov2(end,end))),':',...
                'Color',clrs{2},'LineWidth',2);
    nfpl(2)=plot(rhovec,normpdf(rhovec,th1(:,3),sqrt(cov3(end,end))),'--',...
                'Color',clrs{3},'LineWidth',2);
    xlabel('$\rho$')

    % Calculations for error ellipses
    t=linspace(0,2*pi);
    cl=0.95;
    s=chi2inv(cl,2);
    [V,D]=eig(matslice(covemp,ifinv));
    aemp=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    [V,D]=eig(cov1);
    asample=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    [V,D]=eig(cov2);
    a2=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    [V,D]=eig(cov3);
    a3=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    % Cross plot with error ellipses
    subplot(223)
    spl(1)=plot(thhats(:,1),thhats(:,3),'ko');
    hold on
    er(1)=plot(aemp(1,:)+mobs(1),aemp(2,:)+mobs(3),'Color','k');
    er(2)=plot(asample(1,:)+mobs(1),asample(2,:)+mobs(3),'Color',clrs{1});
    er(3)=plot(a2(1,:)+mobs(1),a2(2,:)+mobs(3),':','Color',clrs{2});
    er(4)=plot(a3(1,:)+mobs(1),a3(2,:)+mobs(3),'--','Color',clrs{3});
    vpl(3)=xline(th1(:,1),'Label','\sigma^2_0','LabelOrientation','horizontal');
    hpl(3)=yline(th1(:,3),'Label','\rho_0','LabelOrientation','horizontal');
    xlabel('$\hat{\sigma^2}$','Interpreter','latex')
    ylabel('$\hat{\rho}$','Interpreter','latex')
    grid on
    box on

    % Covariance comparison
    subplot(224)
    spl(2)=plot(1:numsims,cov1s2s2,'-','Color',clrs{1},'DisplayName','{s2,s2}');
    hold on
    spl(3)=plot(1:numsims,cov1s2rh,'--','Color',clrs{1},'DisplayName','{s2,rh}');
    spl(4)=plot(1:numsims,cov1rhrh,':','Color',clrs{1},'DisplayName','{rh,rh}');
    hJ(1)=yline(covemp(1,1),'-','Color','k');
    hJ(2)=yline(covemp(1,end),'--','Color','k');
    hJ(3)=yline(covemp(end,end),':','Color','k');
    hJ(4)=yline(cov2(1,1),'-','Color',clrs{2});
    hJ(5)=yline(cov2(1,end),'--','Color',clrs{2});
    hJ(6)=yline(cov2(end,end),':','Color',clrs{2});
    hJ(7)=yline(cov3(1,1),'-','Color',clrs{3});
    hJ(8)=yline(cov3(1,end),'--','Color',clrs{3});
    hJ(9)=yline(cov3(end,end),':','Color',clrs{3});
    hJ(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(3).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(4).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(5).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(6).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(7).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(8).Annotation.LegendInformation.IconDisplayStyle = 'off';
    hJ(9).Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlabel('number of simulations')
    ylabel('cov')
    grid on;
    legend()

    sgtitle(sprintf('Exponential Model, \\sigma^2=%0.2f, \\rho=%0.2f; %im x %im',th1(1),th1(3),params.NyNx.*params.dydx))
    saveas(gcf,...
      sprintf('covgammiosl_demo2_2pcrossplots_%ix%i_s2%inu%irh%i_%s.eps',...
        params.NyNx.*params.dydx,th0.*[1 10 1],date),'epsc')

    % Normalized covariance heatmaps
    figure()
    kelicol
    cmrange=[0.5 1];
    labs={'s2','nu','rho'};
    labs=matslice(labs,ifinv);
    np=sum(ifinv);
    ah(1)=subplot(221);
    im(1)=imagesc(matslice(covemp,ifinv)./...
          [diag(sqrt(matslice(covemp,ifinv)))*...
           diag(sqrt(matslice(covemp,ifinv)))'],...
           cmrange);
    title('empirical')
    axis image

    ah(2)=subplot(222);
    im(2)=imagesc(cov1./[diag(sqrt(cov1))*diag(sqrt(cov1))'],cmrange);
    title('sample')
    axis image

    ah(3)=subplot(223);
    im(3)=imagesc(cov2./[diag(sqrt(cov2))*diag(sqrt(cov2))'],cmrange);
    title('diagonals')
    axis image

    ah(4)=subplot(224);
    im(4)=imagesc(cov3./[diag(sqrt(cov3))*diag(sqrt(cov3))'],cmrange);
    title('dftmtx')
    axis image

    set(ah(:),'xtick',1:np,'XTickLabel',labs)
    set(ah(:),'ytick',1:np,'YTickLabel',labs)
    longticks(ah(:))

    cb=colorbar;
    cb.Label.String='normalized';

    sgtitle(sprintf('Exponential Model, \\sigma^2=%0.2f, \\rho=%0.2f; %im x %im',th1(1),th1(3),params.NyNx.*params.dydx))
    saveas(gcf,...
      sprintf('covgammiosl_demo2_2pcorrplots_%ix%i_s2%inu%irh%i_%s.eps',...
        params.NyNx.*params.dydx,th0.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo3')
    % Demonstration of SAMPLE method approaching DFTMTX and
    % DIAGONALS methods by increasing the number of periodograms 

    % Option to either have a quick look (varofcov=0) or repeat calculations N
    % times to find the variance of the covariance (varofcov=1)
    varofcov=1;N=100;

    % Set the Matern and grid parameters
    th1=[5 0.5 5]; ifinv=[1 1 1]; np=sum(ifinv);
    params=[];params.NyNx=[58 53];params.dydx=[1 1];params.blurs=-1;params.taper=1; 

    numsims=1e3;
    if varofcov
      cov1=zeros(6,numsims,N);
      for jnd=1:N
        % Output the gradients calculated by the SAMPLES method
        [covg1,grads]=covgammiosl(th1,params,1,ifinv);
        % Calculate the Fisher inverse
        k=knums(params);
        F1=fishiosl(k,th1,params,0);
        F1=matslice(F1,ifinv);
        covF=inv(F1);
        % The inside of COVTHOSL
        for ind=2:numsims
            covg1=real(cov(grads(1:ind,:)))./prod(params.NyNx).^2;
            covg1=matslice(covg1,[ifinv]);
            cov1(:,ind,jnd)=trilos(covF*covg1*covF)';
        end
      end
      % Calculate the variance for each parameter covariance (3rd dimension)
      [var1,mn1]=var(cov1,0,3);
    else
      % Output the gradients calculated by the SAMPLES method
      [covg1,grads]=covgammiosl(th1,params,1,ifinv);
      % Calculate the Fisher inverse
      k=knums(params);
      F1=fishiosl(k,th1,params,0);
      F1=matslice(F1,ifinv);
      covF=inv(F1);
      cov1=zeros(6,numsims);
      % The inside of COVTHOSL
      for ind=2:numsims
          covg1=real(cov(grads(1:ind,:)))./prod(params.NyNx).^2;
          covg1=matslice(covg1,[ifinv]);
          cov1(:,ind)=trilos(covF*covg1*covF)';
      end
    end

    covg2=covgammiosl(th1,params,2,ifinv);
    covg3=covgammiosl(th1,params,3,ifinv);
    cov2=covthosl(th1,params,covg2,ifinv);
    cov2=trilos(cov2);
    cov3=covthosl(th1,params,covg3,ifinv);
    cov3=trilos(cov3);

    % Labels
    labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
          '\sigma^2,\rho', '\nu,\rho'};
    % Color friendly palette
    clrs=["#1B9E77","#D95F02","#7570B3"];

    clf;
    [ah,ha,H]=krijetem(subnum(2,3));
    for ind=1:6
      axes(ah(ind))
      if varofcov
        pl1(ind)=errorbar(2:25:numsims,mn1(ind,2:25:end),var1(ind,2:25:end),...
          "o",'MarkerSize',2,'Color',clrs{1},'DisplayName','sample');
        pl1(ind).CapSize = 0;
        set(gca,'YScale','log')
      else
        pl1(ind)=plot(2:numsims,cov1(ind,2:end),'-','Color',clrs{1},...
          'LineWidth',2,'DisplayName','sample');
      end
      hold on
      % pl2(ind)=yline(covemp(ind),'-','Color','k');
      pl3(ind)=yline(cov2(ind),':','Color',clrs{2},'LineWidth',2,...
                 'DisplayName','dftmtx');
      pl4(ind)=yline(cov3(ind),'--','Color',clrs{3},'LineWidth',2,...
                 'DisplayName','diagonals');
      if ind==1 | ind==4
        ylabel('Covariance');
      end
      if ind>3
        xlabel('Number of simulations');
      end
      grid on; longticks 
      if ind==1
        lg1=legend();
        lg1.AutoUpdate='off';
      end
      title(sprintf('%s',labs{ind}))
      voff=max(10^floor(log10(abs(cov2(ind)))),max(var1(ind,:))/2);
      ylim([cov2(ind)-voff cov2(ind)+voff])
    end

   keyboard
   ti=sgtitle(sprintf('\\sigma^2=%0.2f, \\nu=%0.2f, \\rho=%0.2f; %im x %im',...
           th1,params.NyNx.*params.dydx))
   movev(ah,-0.03)

   saveas(gcf,...
     sprintf('covgammiosl_demo3_sample_%ix%i_s2%inu%irh%i_%s.eps',...
       params.NyNx.*params.dydx,th1.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo4')
    % Visualization of the construction of terms for the DFTMTX method (distance 
    % grid, covariance, Hk2cov) for very (unrealistically) small grids
  
    params=[];
    params.NyNx=[9 9];
    params.dydx=[1 1];
    params.blurs=-1;
    params.taper=1;

    ifinv=[1 1 1];
    th=[1 0.5 2];
  
    dydx=params.dydx;NyNx=params.NyNx;
    ys=0:NyNx(1)-1;xs=0:NyNx(2)-1;
    [grdx,grdy]=meshgrid(xs(:),ys(:));
    grdx=grdx(:).*dydx(2);
    grdy=grdy(:).*dydx(1);
    lagx=grdx-grdx';
    lagy=grdy-grdy';
    dxxp=sqrt(lagx.^2+lagy.^2);
    Cmn=maternosy(dxxp,th);
    Cm_n_mn=reshape(Cmn,[NyNx(1) NyNx(2) prod(NyNx)]);
    Um=dftmtx(size(Cm_n_mn,1));
    Un=dftmtx(size(Cm_n_mn,2));
    EJJht_in=reshape(permute(...
        tensorprod(Um,tensorprod(Un,Cm_n_mn,1,2),1,2),...
        [3 1 2]),NyNx(1),NyNx(2),prod(NyNx));
    % This will be going to the "Sbar"
    EJJht_out=reshape(permute(...
                conj(tensorprod(Um,tensorprod(Un,conj(EJJht_in),1,2),1,2)),...
              [3 1 2]),[prod(NyNx) prod(NyNx)]);
    % This will be going to the "pSbar"
    EJJt_out=reshape(permute(...
               tensorprod(Um,tensorprod(Un,EJJht_in,1,2),1,2),...
             [3 1 2]),[prod(NyNx) prod(NyNx)]);
    Hk2cov=(abs(EJJht_out).^2+abs(EJJt_out).^2)./prod(NyNx).^2;
  
    % While the first and second terms below sum to similar quantities, we can
    % quickly see wave-number dependent distinctions:
    figure(1)
    clf 
    ah(1)=subplot(321); imagesc(dxxp); cb1=colorbar; 
    cb1.Label.String='euclidean distance [m]'; 
    t(1)=title('$||\mathbf{x}-\mathbf{x''}||$');
    t(1).Interpreter='latex';
    % Ticks represent indices only
    ah(1).XTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    ah(1).YTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    xlabel('$\mathbf{x}''$'); ylabel('$\mathbf{x}$'); axis image

    ah(2)=subplot(322); imagesc(Cmn); caxis([0 1]); cb2=colorbar; 
    cb2.Label.String='covariance [unit^2]'; 
    t(2)=title('$C(||\mathbf{x}-\mathbf{x''}||)$');
    t(2).Interpreter='latex';
    ah(2).XTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    ah(2).YTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    xlabel('$\mathbf{x}''$'); ylabel('$\mathbf{x}$'); axis image

    ah(3)=subplot(323);
    imagesc(log(abs(EJJt_out./prod(NyNx)).^2));
    cax=caxis; caxis(cax(2)-[10 0])
    cb4=colorbar; cb4.Label.String='term 1 log([|unit|^2/m])';
    %title('|U(UCU'')U|^2')
    t(3)=title('$|\mathrm{cov}\{H(\mathbf{k}),H^*(\mathbf{k''})\}|^2$');
    t(3).Interpreter='latex';
    ah(3).XTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    ah(3).YTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    xlabel('$\mathbf{k}''$'); ylabel('$\mathbf{k}$'); axis image

    ah(4)=subplot(324);
    imagesc(log(abs(EJJht_out./prod(NyNx)).^2));
    cax=caxis; caxis(cax(2)-[10 0])
    cb3=colorbar; cb3.Label.String='term 2 log([|unit|^2/m])';
    %title('|(U(UCU'')*U)*|^2')
    t(4)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k''})\}|^2$');
    t(4).Interpreter='latex';
    ah(4).XTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    ah(4).YTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    xlabel('$\mathbf{k}''$'); ylabel('$\mathbf{k}$'); axis image

    % Ratio of Isserlis' terms for periodogram covariance
    ah(5)=subplot(325);
    imagesc(log(Hk2cov));
    cax=caxis; caxis(cax(2)-[10 0])
    cb5=colorbar; cb5.Label.String='log periodogram covariance';
    t(5)=title('$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k''})|^2\}$');
    t(5).Interpreter='latex';
    ah(5).XTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    ah(5).YTick=[0:params.NyNx(1):prod(params.NyNx)-1 prod(params.NyNx)];
    xlabel('$\mathbf{k}''$'); ylabel('$\mathbf{k}$'); axis image
    layout(ah(5),0.5,'m','x')
    
    % Check this out this should be the covariance between periodogram terms
    A=Hk2cov;

    % %clf
    % % Grab a 3x3 piece using BLOCK
    % %imagesc(log10(fftshift(blocktile(A,9,0,randi(81)))))
    % indi=0;
    % % step down the block diagonal
    % for index=1+[NyNx(1)+1]*[0:NyNx(1)-1]
    %     indi=indi+1;
    %     % Attempt to find the k=k'
    %     mm(indi)=max(max(fftshift(blocktile(A,NyNx(1),0,index))));
    % end

    % The first term of this is equal to mm(1), but later terms are larger
    % because we are not accounting for the anti-diagonal structure of EJJt_out
    Sbar=abs(blurosy(th,params,0,'efd',[0 0]))*(2*pi)^2;
    Sbar2=abs(Sbar).^2;

    % This is a subset of the k=kp but only one of them
    %mm.'./(diag(v2s(Sbar2))*2)
    % Similarly, the diagonal of A is only twice Sbar^2 for the 0,0 element
    %diag(A)./(Sbar2*2)

    % Hk2cov is the sum of two periodogram terms. The first, abs(EJJht_out).^2,
    % is equivalent to the periodogram calculated from the blurred, shifted 
    % spectral density with no offset
    difer(diag(EJJht_out)./prod(NyNx)-Sbar)

    % The second periodogram term we need to predict is EJJt_out
    pSbar=zeros(prod(NyNx),1);
    % Assumes square grid, Ny==Nx
    if mod(NyNx(1),2)
      % odd
      seq=[0:2:floor(NyNx(1)/2).*2 -floor(NyNx(1)/2).*2:2:-2];
    else
      % even
      seq=[0:2:NyNx(1) -NyNx(1)+2:2:-2];
    end
    % this running index...
    rnd=0;
    for jnd=seq
      for ind=seq
        rnd=rnd+1;
        Sbaroff=blurosy(th,params,0,'efd',[ind jnd])*(2*pi)^2;
        % ... finds the relevant k=k' term
        pSbar(rnd,:)=Sbaroff(rnd);
      end
    end
    difer(diag(EJJt_out)./prod(NyNx)-pSbar)

    % The diagonal of the sum of the two predicted terms is the diagonal of A
    pSbar2=abs(pSbar).^2;

    diagpgmsum=v2s(Sbar2+pSbar2);
    difer(diag(A)-diagpgmsum(:))
    % Bottom line: Hk2cov is the sum of a covariance of Fourier coefficients and
    % a pseudocovariance of Fourier coefficients. The variances, aka their
    % diagonals, can be predicted from the blurred spectral density.

    figure(2)
    clf
    
    % This is more than we wanted to know but might some some day
    % The first column (or equivalently the first row) of the wrapped diagonal
    % of the sum of the predicted terms is mm
    % difer(mm.'-diagpgmsum(:,1))
    % Don't v2s(fftshift, very different from fftshitf(v2s
    ah(1)=subplot(221); imagesc(decibel(fftshift(v2s(pSbar2))))
    caxis([-40 0])
    t(1)=title('$|\mathrm{cov}\{H(\mathbf{k}),H^*(\mathbf{k})\}|^2$');
    t(1).Interpreter='latex';
    ah(1).XTick=[1 params.NyNx(1)];
    ah(1).YTick=[1 params.NyNx(2)];
    xlabel('$k_x$'); ylabel('$k_y$'); axis image

    ah(2)=subplot(222); imagesc(decibel(fftshift(v2s(Sbar2))))
    caxis([-40 0])
    t(2)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k})\}|^2$');
    t(2).Interpreter='latex';
    % Ticks represent indices only
    ah(2).XTick=[1 params.NyNx(1)];
    ah(2).YTick=[1 params.NyNx(2)];
    xlabel('$k_x$'); ylabel('$k_y$'); axis image

    a=fftshift(v2s(Sbar2+pSbar2));
    ah(3)=subplot(223); imagesc(decibel(a));
    caxis([-40 0])
    t(3)=title('$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k''})|^2\}$');
    t(3).Interpreter='latex';
    ah(3).XTick=[1 params.NyNx(1)];
    ah(3).YTick=[1 params.NyNx(2)];
    xlabel('$k_x$'); ylabel('$k_y$'); axis image

    % Note that this is only the non-pseudo part of what we need
    % The non-correlated part of the covariance of the periodogram is not just
    % the square of the blurred theoretical spectrum
    b=v2s(blurosy(th,params)*(2*pi)^2).^2;
    [~,i]=max(b(:)); b(i)=b(i)*2;
    ah(4)=subplot(224); imagesc(decibel(b))
    caxis([-40 0])
    t(4)=title(sprintf('fixed %s [median ratio %5g]','$\bar{S}(\mathbf{k})^2$',median(abs(a(:)./b(:)))));
    t(4).Interpreter='latex';
    ah(4).XTick=[1 params.NyNx(1)];
    ah(4).YTick=[1 params.NyNx(2)];
    xlabel('$k_x$'); ylabel('$k_y$'); axis image

    % layout(ah(3),0.5,'m','x')
    
    % % Check if Hermitian
    % Sbar2=fftshift(v2s(Sbar2));
    % pSbar2=fftshift(v2s(pSbar2));
    % hermcheck(Sbar2)
    % hermcheck(pSbar2)
    % % Check if positive definite 
    % flg=defcheck(Sbar2); diferm(flg>1)
    % pflg=defcheck(pSbar2); diferm(pflg>1)

    keyboard
    saveas(gcf,...
      sprintf('covgammiosl_demo4_dftmtx_%ix%i_s2%inu%irh%i_%s.eps',...
        params.NyNx.*params.dydx,th.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo5')
   % [IN PROGRESS] Visual demonstration of the indexing scheme implemented in the DIAGONALS
   % method
   % TODO: COPY IN SCRATCH/OLIVIASJ.M
elseif strcmp(th,'demo6')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the GRID SIZE, which is presented in terms
   % of pi rho by pi rho. Plot the bias from the empirical study for each
   % parameter.
   params=[];params.dydx=[1 1];params.blurs=-1;params.taper=1;ifinv=[1 1 1];
   numreals=200;

   % Labels
   labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
         '\sigma^2,\rho', '\nu,\rho'};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];

   % For now, small values so that we have a chance at forming the data vector
   % and finding an estimate for small grids 
   rh=2;
   th1=[1 0.5 rh];
   % Sample in logspace
   loggs=logspace(0,1,10);
   thhats=zeros(numreals,3,size(loggs,2));
   try
     % If we have already calculate the estimates, load them. Note that we may
     % need to append a date to the file name
     datum='26-Mar-2025';%date;
     thhats=load(sprintf('covgammiosl_demo6_empthhats_%s.mat',datum),'thhats');
     thhats=thhats.thhats;
   catch
     for ind=1:size(loggs,2)
       % Create a larger grid with length linearly increasing in pi rho
       gs=rh*pi*loggs(ind); 
       params.NyNx=[floor(gs) floor(gs)];
       % Simulate via circulant embedding
       simparams=params;simparams.blurs=Inf;
       for mnd=1:numreals
         Hx=simulosl(th1,simparams,0);
         thhat=NaN;
         while isnan(thhat)
           [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
         end
         thhats(mnd,:,ind)=thhat.*scl;
       end
     end
     keyboard
     for ind=1:size(loggs,2)
       % Consider trimming, some results are truly bad on small grids
       [thhatst{ind},trimi]=trimit(thhats(:,:,ind),97);
       for jnd=1:(numreals-size(trimi,1))
         Hx=simulosl(th1,simparams,0);
         thhat=NaN;
         while isnan(thhat)
           [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
         end
         thhats2(jnd,:)=thhat.*scl;
       end
       % Recreate thhats(:,:,ind) with the recalculated estimates
       thhats(:,:,ind)=[thhatst{ind};thhats2];
     end
     % Consider saving your calculations by setting 'svthhats' in the debugger
     keyboard
     svthhats=1;
     if svthhats
       save(sprintf('covgammiosl_demo6_empthhats_%s.mat',date),'thhats')
     end
   end

   % Form the covariance matrices
   for ind=1:size(loggs,2)
     gs=rh*pi*loggs(ind);
     params.NyNx=[floor(gs) floor(gs)];
     % Calculate the empirical covariance from the ensemble
     covemp=nancov(thhats(:,:,ind));
     % Calculate the 'sample' covariance
     covg1  =covgammiosl(th1,params,1,ifinv);
     cov1c=covthosl(th1,params,covg1,ifinv);
     % Calculate the 'dftmtx' covariance
     covg2  =covgammiosl(th1,params,2,ifinv);
     cov2c=covthosl(th1,params,covg2,ifinv);
     % Calculate the 'diagonals' covariance
     covg3  =covgammiosl(th1,params,3,ifinv);
     cov3c=covthosl(th1,params,covg3,ifinv);
     % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
     covempt(ind,:)=trilos(covemp);
     cov1t(ind,:)  =trilos(cov1c);
     cov2t(ind,:)  =trilos(cov2c);
     cov3t(ind,:)  =trilos(cov3c);
   end
   cmnylim=[minmax([abs(covempt) abs(cov1t) abs(cov2t) abs(cov3t)])];
   cmnylimrng=round(log10(cmnylim),0);
   clf;
   [ah,ha,H]=krijetem(subnum(3,3));
   for jnd=1:6
     axes(ah(jnd))
     % Plot reference diagonals as light gray lines on bottom for each power of
     % 10 between cmnylim
     hold on
     rnd=0;
     for ind=cmnylimrng(1):cmnylimrng(2)+1
       rnd=rnd+1;
       refd(rnd)=loglog(loggs,10^(ind)./floor(rh.*pi.*loggs),...
                 'Color',[0.8 0.8 0.8],'HandleVisibility','off');
     end
     % Fit a slope to the covariance calculations as a function of grid length
     % in log-log
     mbfemp(jnd,:) =polyfit(log10(loggs),log10(abs(covempt(:,jnd))),1);
     mbfcov1(jnd,:)=polyfit(log10(loggs),log10(abs(cov1t(:,jnd))),1);
     mbfcov2(jnd,:)=polyfit(log10(loggs),log10(abs(cov2t(:,jnd))),1);
     mbfcov3(jnd,:)=polyfit(log10(loggs),log10(abs(cov3t(:,jnd))),1);
     predemp(jnd,:) =polyval(mbfemp(jnd,:), log10(loggs));
     predcov1(jnd,:)=polyval(mbfcov1(jnd,:),log10(loggs));
     predcov2(jnd,:)=polyval(mbfcov2(jnd,:),log10(loggs));
     predcov3(jnd,:)=polyval(mbfcov3(jnd,:),log10(loggs));
     % Plot the absolute covariance so that negative correlations do not 
     % mislead
     ll1(ind)=loglog(loggs,abs(covempt(:,jnd)),'ko',...
            'DisplayName',sprintf('empirical, %0.1f',abs(mbfemp(jnd,1))));
     ll2(ind)=loglog(loggs,abs(cov1t(:,jnd)),'^','MarkerEdgeColor',clrs{1},...
            'MarkerFaceColor',clrs{1},...
            'DisplayName',sprintf('sample, %0.1f',abs(mbfcov1(jnd,1))));
     ll3(ind)=loglog(loggs,abs(cov2t(:,jnd)),'s','MarkerEdgeColor',clrs{2},...
            'MarkerFaceColor',clrs{2},...
            'DisplayName',sprintf('dftmtx, %0.1f',abs(mbfcov2(jnd,1))));
     ll4(ind)=loglog(loggs,abs(cov3t(:,jnd)),'*','MarkerEdgeColor',clrs{3},...
            'MarkerFaceColor',clrs{3},...
            'DisplayName',sprintf('diagonals, %0.1f',abs(mbfcov3(jnd,1))));
     % Plot the best fit lines
     loglog(loggs,10.^(predemp(jnd,:)),'Color','k','HandleVisibility','off');
     loglog(loggs,10.^(predcov1(jnd,:)),'Color',clrs{1},'HandleVisibility','off');
     loglog(loggs,10.^(predcov2(jnd,:)),'Color',clrs{2},'HandleVisibility','off');
     loglog(loggs,10.^(predcov3(jnd,:)),'Color',clrs{3},'HandleVisibility','off');
   end
   for jnd=1:6
     axes(ah(jnd));
     % List slopes in the legend and explain in a figure caption
     [leg1(jnd),legic]=legend('interpreter','tex','BackgroundAlpha',0,'box','off');
     % [leg1(jnd),legic]=legend([ll1(jnd) ll2(jnd) ll3(jnd) ll4(jnd)],...
     %                  'interpreter','tex','BackgroundAlpha',0,'box','off');
     % Have the legend markers take up less space
     legicms = findobj(legic,'Type','line');
     for ind=1:size(legicms,1)
       if size(legicms(ind).XData,1)==1
          legicms(ind).XData=0.3;
       else
          legicms(ind).XData(1)=0.2;
          legicms(ind).XData(2)=0.4;
       end
     end
     % Resize the marker in the legend
     
     % Labels, titles, ticks, common limits
     if jnd==1|jnd==4
        ylabel('$|\mathrm{cov}\{\theta_i,\theta_j\}|$','Interpreter','latex')
     end
     if jnd>3
     %   xlabel('Length ($\pi\rho$)','Interpreter','latex')
     end
     title(sprintf('%s',labs{jnd}))
     longticks
     ah(jnd).XRuler.TickLabelGapOffset=-2;
     box on
     ylim([cmnylim(1) 100*cmnylim(2)])
     xlim([1 10])
     set(gca,'Yscale','log','Xscale','log')
   end      
    
   % Plot the bias of the estimates from the empirical study
   axes(ah(7));
   % semilogx(loggs,squeeze(mean(thhats(:,1,:))-th1(1)),'ko')
   % Errorbars in terms of standard error for the number of realizations
   pl1(1)=errorbar(loggs,squeeze(mean(thhats(:,1,:)-th1(1))),...
                   squeeze(sqrt(var(thhats(:,1,:)-th1(1))/numreals)),'ko');
   pl1(1).CapSize=0;
   % Add the 0-line where the estimate equals the truth
   yline(0)
   ylabel('$\langle\hat{\theta}\rangle-\theta_0$','interpreter','latex');
   ti(7)=title('\sigma^2','FontWeight','bold','Interpreter','tex');
   %
   axes(ah(8));
   % semilogx(loggs,squeeze(mean(thhats(:,2,:))-th1(2)),'ko')
   pl1(2)=errorbar(loggs,squeeze(mean(thhats(:,2,:)-th1(2))),...
                   squeeze(sqrt(var(thhats(:,2,:)-th1(2))/numreals)),'ko');
   pl1(2).CapSize=0;
   yline(0)
   ti(8)=title('\nu','FontWeight','bold','Interpreter','tex');
   ylim([-0.1,max(mean(thhats(:,2,:))-th1(2))+0.1])
   %
   axes(ah(9));
   pl1(3)=errorbar(loggs,squeeze(mean(thhats(:,3,:)-th1(3))),...
                   squeeze(sqrt(var(thhats(:,3,:)-th1(3))/numreals)),'ko');
   pl1(3).CapSize=0;
   yline(0)
   ti(9)=title('\rho','FontWeight','bold','Interpreter','tex');
   %
   for jnd=7:9
     axes(ah(jnd))
     longticks
     xlabel('Length ($\pi\rho$)','Interpreter','latex')
     ah(jnd).XRuler.TickLabelGapOffset=-2;
     set(gca,'XScale','log')
   end
   %

   % Super title to display the parameters
   ti=sgtitle(...
              sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s','\sigma^2=',th1(1),...
                  ' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),' m'),...
              'Interpreter','latex');
   % Shift the figure objects to make space for the title
   movev([ah(1:3) leg1(1) leg1(2) leg1(3)],-0.03)
   movev([ah(4:6) leg1(4) leg1(5) leg1(6)],-0.01)
   movev(ah(7:9),0.01)

   % Save the figure
   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo6_gridsize_s2%inu%irh%i_%s.eps',...
       th1.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo7')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the GRID SPACING. This is a numerical
   % implementation of studying INFILL ASYMPTOTICS (Cressie, 1993; also referred
   % to as FIXED-DOMAIN ASYMPTOTICS by Stein, 1999), where sample size is
   % increased within a bounded domain (i.e., NyNx.*dydx is fixed)
   params=[];params.blurs=-1;params.taper=1;ifinv=[1 1 1];
   numreals=100;

   % Labels
   labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
          '\sigma^2,\rho', '\nu,\rho'};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];

   % Set the Matern parameters 
   rh=2;
   th1=[1 0.5 rh];
   % Make a square grid with fixed grid length, gy
   gy=120;
   gygx=[gy gy];
   % Calculate factors of gy to choose multiple instances of NyNx and dydx with
   % the same gygx
   facs=1:ceil(sqrt(gy));
   dy=facs(rem(gy,facs)==0);
   ny=sort(gy./dy,'descend'); 
   % Confirm that each set of ny and dy results in the fixed grid size
   diferm(ny.*dy-gy)
   % We should only take stock from experiments that have a sufficiently large
   % number of samples to begin with, which we will impose as 3 times pi*rho
   % (still pretty small)
   dy=dy(ny>=pi*rh*3);
   ny=ny(ny>=pi*rh*3);
   % This looks the nicest for gy==120
   if gy==120
     dy=[1 2 3 4 6];
     ny=[120 60 40 30 20];
   end

   % We need an MLEOSL suite
   try
     % In case we already calculated the estimates, we will try loading them
     datum=date;
     thhats=load(sprintf('covgammiosl_demo7_empthhats_%s.mat',datum),'thhats');
     thhats=thhats.thhats;
     if gy==120&strcmp(datum,'07-Apr-2025')
       thhats=thhats(:,:,[1:4,6]);
     end
     keyboard
   catch
     for ind=1:size(dy,2)
       % Set spacing and number of samples
       params.NyNx=[ny(ind) ny(ind)];
       params.dydx=[dy(ind) dy(ind)];
       % Simulate via circulant embedding
       simparams=params;simparams.blurs=Inf;
       for mnd=1:numreals
         Hx=simulosl(th1,simparams,0);
         thhat=NaN;
         while isnan(thhat)
           [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
         end
         thhats(mnd,:,ind)=thhat.*scl;
       end
     end
     % Consider saving your calculations by setting 'svthhats' in the debugger
     keyboard
     svthhats=1;
     if svthhats
       save(sprintf('covgammiosl_demo7_empthhats_%s.mat',date),'thhats')
     end
   end
   % Form the covariance matrices
   for ind=1:size(dy,2)
     % Set spacing and number of samples
     params.NyNx=[ny(ind) ny(ind)];
     params.dydx=[dy(ind) dy(ind)];
     % Calculate the empirical covariance from the ensemble of estimates
     covemp=nancov(thhats(:,:,ind));
     % Calculate parameter covariance according to the 'sample' method
     covg1  =covgammiosl(th1,params,1,ifinv);
     cov1c=covthosl(th1,params,covg1,ifinv);
     % Calculate parameter covariance according to the 'dftmtx' method
     covg2  =covgammiosl(th1,params,2,ifinv);
     cov2c=covthosl(th1,params,covg2,ifinv);
     % Calculate parameter covariance according to the 'diagonals' method
     covg3  =covgammiosl(th1,params,3,ifinv);
     cov3c=covthosl(th1,params,covg3,ifinv);
     % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
     covempt(ind,:)=trilos(covemp);
     cov1t(ind,:)  =trilos(cov1c);
     cov2t(ind,:)  =trilos(cov2c);
     cov3t(ind,:)  =trilos(cov3c);
   end

   clf;
   [ah,ha,H]=krijetem(subnum(3,3));
   cmnylim=[minmax([abs(covempt) abs(cov1t) abs(cov2t) abs(cov3t)])];
   cmnylimrng=ceil(log10(cmnylim));
   for jnd=1:6
     axes(ah(jnd))
     % Plot reference diagonals as light gray lines on bottom for each power of
     % 10 between cmnylim
     hold on
     rnd=0;
     for ind=cmnylimrng(1):cmnylimrng(2)
       rnd=rnd+1;
       refd(rnd)=loglog(ny,10^(ind)./ny,'Color',[0.8 0.8 0.8]);
     end
     % Plot the absolute covariance so that negative correlations do not 
     % mislead
     loglog(ny,abs(covempt(:,jnd)),'ko','DisplayName','empirical')
     loglog(ny,abs(cov1t(:,jnd))  ,'^','MarkerEdgeColor',clrs{1},...
            'MarkerFaceColor',clrs{1},'DisplayName','sample')
     loglog(ny,abs(cov2t(:,jnd))  ,'s','MarkerEdgeColor',clrs{2},...
            'MarkerFaceColor',clrs{2},'DisplayName','dftmtx')
     loglog(ny,abs(cov3t(:,jnd))  ,'*','MarkerEdgeColor',clrs{3},...
            'MarkerFaceColor',clrs{3},'DisplayName','diagonals')
   end

   for jnd=1:6
     axes(ah(jnd));
     % Fit a slope to the covariance calculations as a function of grid length
     % in log-log
     mbfemp(jnd,:) =polyfit(log10(ny),log10(abs(covempt(:,jnd))),1);
     mbfcov1(jnd,:)=polyfit(log10(ny),log10(abs(cov1t(:,jnd))),1);
     mbfcov2(jnd,:)=polyfit(log10(ny),log10(abs(cov2t(:,jnd))),1);
     mbfcov3(jnd,:)=polyfit(log10(ny),log10(abs(cov3t(:,jnd))),1);
     predemp(jnd,:) =polyval(mbfemp(jnd,:), log10(ny));
     predcov1(jnd,:)=polyval(mbfcov1(jnd,:),log10(ny));
     predcov2(jnd,:)=polyval(mbfcov2(jnd,:),log10(ny));
     predcov3(jnd,:)=polyval(mbfcov3(jnd,:),log10(ny));
     loglog(ny,10.^(predemp(jnd,:)),'Color','k');
     loglog(ny,10.^(predcov1(jnd,:)),'Color',clrs{1});
     loglog(ny,10.^(predcov2(jnd,:)),'Color',clrs{2});
     loglog(ny,10.^(predcov3(jnd,:)),'Color',clrs{3});

     [leg1(jnd),legic]=legend([repmat("",rnd,1)',...
                      sprintf('empirical, %0.1f',abs(mbfemp(jnd,1))),...
                      sprintf('sample, %0.1f',abs(mbfcov1(jnd,1))),...
                      sprintf('dftmtx, %0.1f',abs(mbfcov2(jnd,1))),...
                      sprintf('diagonal, %0.1f',abs(mbfcov3(jnd,1)))],...
                      'interpreter','tex','BackgroundAlpha',0,'box','off');
     leg1(jnd).AutoUpdate='off';
     % Have the legend markers take up less space
     legicms = findobj(legic,'Type','line');
     for ind=1:size(legicms,1)
       if size(legicms(ind).XData,1)==1
          legicms(ind).XData=0.3;
       else
          legicms(ind).XData(1)=0.2;
          legicms(ind).XData(2)=0.4;
       end
     end
     % Resize the marker in the legend
     
     % Labels, titles, ticks, common limits
     if jnd==1|jnd==4
        ylabel('$|\mathrm{cov}\{\theta_i,\theta_j\}|$','Interpreter','latex')
     end
     if jnd>3
     %   xlabel('Length ($\pi\rho$)','Interpreter','latex')
     end
     ti(jnd)=title(sprintf('%s',labs{jnd}));
     longticks
     ylim([cmnylim(1)*0.5 cmnylim(2)*2.0])
     longticks
     box on
     xlim(minmax(ny).*[0.9 1.1])
     set(gca,'Yscale','log','Xscale','log','XTick',fliplr(ny),'XTickLabel',...
         fliplr(ny))
     % We will want a second x-axis that quotes the spacing (dy); adapted from
     % eggers2
     xlb='Spacing (\Deltax)';
     if jnd<4
       [axx(jnd),xl(jnd),yl(jnd)]=...
         xtraxis(ah(jnd),dy,dy,xlb);
     %  movev(ti(jnd),15);
     else
       [axx(jnd),xl(jnd),yl(jnd)]=...
         xtraxis(ah(jnd),dy,dy,[]);
     end
     axx(jnd).XLim=fliplr(gy./minmax(ny).*[1.1 0.9]);
     set(axx(jnd),'xdir','rev')
     longticks
     %movev(xl(jnd),-10)
     ah(jnd).XRuler.TickLabelGapOffset=-2;
     ah(jnd).XTickLabelRotation=0;
     if gy==120
       ah(jnd).XTickLabel(3,:)=' ';
     end
     axx(jnd).XRuler.TickLabelGapOffset=-3;
   end      
    
   xla='Number of Samples';
   % Plot the bias of the estimates from the empirical study
   axes(ah(7));
   % Errorbars in terms of standard error for the number of realizations
   pl1(1)=errorbar(ny,squeeze(mean(thhats(:,1,:)-th1(1))),...
                   squeeze(sqrt(var(thhats(:,1,:)-th1(1))/numreals)),'ko');
   pl1(1).CapSize = 0;
   % Add the 0-line where the estimate equals the truth
   yline(0)
   xlabel(xla)
   ylabel('$\langle\hat{\theta}\rangle-\theta_0$','interpreter','latex');
   ti(7)=title('\sigma^2','FontWeight','bold','Interpreter','tex');
   %
   axes(ah(8));
   pl1(2)=errorbar(ny,squeeze(mean(thhats(:,2,:)-th1(2))),...
                   squeeze(sqrt(var(thhats(:,2,:)-th1(2))/numreals)),'ko');
   pl1(2).CapSize = 0;
   yline(0)
   xlabel(xla)
   ti(8)=title('\nu','FontWeight','bold','Interpreter','tex');
   ylim([-0.1,max(mean(thhats(:,2,:))-th1(2))+0.1])
   %
   axes(ah(9));
   pl1(3)=errorbar(ny,squeeze(mean(thhats(:,3,:)-th1(3))),...
                   squeeze(sqrt(var(thhats(:,3,:)-th1(3))/numreals)),'ko');
   pl1(3).CapSize = 0;
   yline(0)
   xlabel(xla)
   ti(9)=title('\rho','FontWeight','bold','Interpreter','tex');
   %
   byl=minmax([ylim(ah(7)) ylim(ah(8)) ylim(ah(9))]);
   for jnd=7:9
     axes(ah(jnd))
     longticks
     xlim(minmax(ny).*[0.9 1.1])
     set(gca,'Xscale','log','XTick',fliplr(ny),'XTickLabel',...  
         fliplr(ny),'YLim',byl)          
     [axx(jnd),xl(jnd),yl(jnd)]=...
       xtraxis(ah(jnd),dy,dy,[]);
     axx(jnd).XLim=fliplr(gy./minmax(ny).*[1.1 0.9]);
     set(axx(jnd),'xdir','rev')
     longticks
     %movev(xl(jnd),-10)
     ah(jnd).XRuler.TickLabelGapOffset=-2;
     ah(jnd).XTickLabelRotation=0;
     axx(jnd).XRuler.TickLabelGapOffset=-3;
   end
   %

   % Super title to display the parameters
   sti=sgtitle(...
              sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s','\sigma^2=',th1(1),...
                  ' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),' m'),...
              'Interpreter','latex');
   % Shift the figure objects to make space for the title
   shrink([ah axx],1.05,1.05)
   movev([axx(1:3) ah(1:3)],-0.07)
   movev([axx(4:6) ah(4:6)],-0.06)
   movev([axx(7:9) ah(7:9)],-0.05)
   % Move the axes titles
   xlpv=xl(1).Position(2);
   tipv=ti(1).Position(2);
   for jnd=1:3
     xl(jnd).Position(2)=5*xlpv;
     ti(jnd).Position(2)=0.75*xlpv;
   end
   for jnd=4:6
     ti(jnd).Position(2)=0.75*xlpv;
   end
   movev([ti(7:9)],max(byl)/10)

   % Save the figure
   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo7_infill_%inu%irh%i_%s.eps',...
       th1.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo8')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the PERCENT OF MISSINGNESS, either due to
   % random deletions or a growing mask

   % The Matern and regular grid parameters
   th1=[3 0.5 4];
   params.NyNx=[88 88];params.dydx=[1 1];
   params.blurs=-1;params.taper=1;                                            
   ifinv=[1 1 1];                                                              
   oneprcnt=ceil(0.01.*prod(params.NyNx));
   % If you specify a mask, it needs to be usable by MASKIT, either as an index
   % matrix or a region name that we have worked with previously
   % params.mask='france';

   % If you do not specify a mask, we will proceed with random deletions
   
   % Labels
   labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
         '\sigma^2,\rho', '\nu,\rho'};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];
   
   clf;
   [ah,ha,H]=krijetem(subnum(3,3));
   % Simulate via circulant embedding
   taper=ones(params.NyNx);
   simparams=params;simparams.blurs=Inf;
   % The number of parameters
   np=3;
   % How many realizations for calculating the empirical covariance?
   numreals=200;
   % What is the maximum percent missingness we will study?
   mostmiss=round(10)+1;
   % Preallocate space for the estimates
   thhats=zeros(numreals,np,mostmiss);
   % Calculate the estimates for increasing missingness for the mask type
   % selected 
   for ind=1:mostmiss
     params.taper=taper;simparams.taper=taper;
     tapersv(:,:,ind)=taper;
     % Calculate the percentage of the grid that is missing
     prcmiss(ind)=1-sum(taper,"all")./prod(params.NyNx);
     for mnd=1:numreals
       Hx=simulosl(th1,simparams,0);
       thhat=NaN;
       while isnan(thhat)
         [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
       end
       thhats(mnd,:,ind)=deal(thhat.*scl);
     end
     % Apply TRIMIT to remove possible errant estimates 
     [thhatst{ind},trimi]=trimit(thhats(:,:,ind),99);
     if isfield(params,'mask')
       [~,I]=maskit(rand(params.NyNx),params,ind/100);
       taper=~I;
     else
       taper(randi(prod(params.NyNx),1,oneprcnt))=deal(0);
     end
   end
   if isfield(params,'mask')
     save(sprintf('covgammiosl_demo8_%s_%s.mat',params.mask,date),...
          'thhatst','tapersv')
   else
     save(sprintf('covgammiosl_demo8_random_%s.mat',date),...
          'thhatst','tapersv')
   end
   % Calculate parameter covariance 4 different ways and plot as a
   % function of the percentage of the regular grid that is missing
   for ind=1:mostmiss
     params.taper=tapersv(:,:,ind);
     covemp=nancov(thhatst{ind});
     covg1 =covgammiosl(th1,params,1,ifinv);
     covg2 =covgammiosl(th1,params,2,ifinv);
     covg3 =covgammiosl(th1,params,3,ifinv);
     cov1c=covthosl(th1,params,covg1,ifinv);
     cov2c=covthosl(th1,params,covg2,ifinv);
     cov3c=covthosl(th1,params,covg3,ifinv);
     cov1t=trilos(cov1c);
     cov2t=trilos(cov2c);
     cov3t=trilos(cov3c);
     covet=trilos(covemp);
     for jnd=1:6
       axes(ah(jnd))
       plot(prcmiss(ind)*100,abs(covet(jnd)),'ko','DisplayName','empirical')
       hold on
       plot(prcmiss(ind)*100,abs(cov1t(jnd)),'^','MarkerEdgeColor',clrs(1),...
            'MarkerFaceColor',clrs(1),'DisplayName','sample')
       plot(prcmiss(ind)*100,abs(cov2t(jnd)),'s','MarkerEdgeColor',clrs(2),...
            'MarkerFaceColor',clrs(2),'DisplayName','dftmtx')
       plot(prcmiss(ind)*100,abs(cov3t(jnd)),'*','MarkerEdgeColor',clrs(3),...
            'MarkerFaceColor',clrs(3),'DisplayName','diagonals')
       if jnd==2 & ind==1
         % Make a legend the first time we touch the top center plot
         [leg1,legic]=legend('empirical','sample','dftmtx','diagonals',...
                             'box','off');
         leg1.AutoUpdate='off';
         % Reduce the amount of space each legend marker uses 
         legicms = findobj(legic,'Type','line');
         for knd=1:size(legicms,1)
           if size(legicms(knd).XData,1)==1
              legicms(knd).XData=0.3;
           else
              legicms(knd).XData(1)=0.2;
              legicms(knd).XData(2)=0.4;
           end
         end
         % Move the legend to the left edge of the axis
         leg1.Position(1)=0.38;
         leg1.Position(2)=0.78;
      end
    end
   end
   for jnd=1:6
    axes(ah(jnd))
    if jnd==1|jnd==4
       ylabel('$|\mathrm{cov}\{\theta_i,\theta_j\}|$','Interpreter','latex')
    end
    set(groot, 'DefaultTextInterpreter', 'latex')
    title(sprintf('$%s$',labs{jnd}))
    longticks
   end      

   % Plot the bias of the estimates from the empirical study
   axes(ah(7));
   % Errorbar in terms of the standard error given the number of realizations
   pl1(1)=errorbar(prcmiss*100,squeeze(mean(thhats(:,1,:)-th1(1))),...
                   squeeze(sqrt(var(thhats(:,1,:)-th1(1))/numreals)),'ko');
   pl1(1).CapSize = 0;
   hold on
   % Add the 0-line where the estimate equals the truth
   yline(0)
   xlabel('Percent missing')
   ylabel('$\langle\hat{\theta}\rangle-\theta_0$','interpreter','latex');
   ti(7)=title('\sigma^2','Interpreter','tex');
   longticks
   %
   axes(ah(8));
   pl1(2)=errorbar(prcmiss*100,squeeze(mean(thhats(:,2,:)-th1(2))),...
                   squeeze(sqrt(var(thhats(:,2,:)-th1(2))/numreals)),'ko');
   pl1(2).CapSize = 0;
   hold on
   yline(0)
   xlabel('Percent missing')
   ti(8)=title('\nu','Interpreter','tex');
   %ylim([-0.1,max(mean(thhats(:,2,:))-th1(2))+0.1])
   longticks
   %
   axes(ah(9));
   pl1(3)=errorbar(prcmiss*100,squeeze(mean(thhats(:,3,:)-th1(3))),...
                   squeeze(sqrt(var(thhats(:,3,:)-th1(3))/numreals)),'ko');
   pl1(3).CapSize = 0;
   hold on
   yline(0)
   xlabel('Percent missing')
   ti(9)=title('\rho','Interpreter','tex');
   longticks
   %

   if isfield(params,'mask')
     masknm=params.mask;
   else
     masknm='random';
   end

   sti(1)=sgtitle(...
     sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s%i%s%i%s%s',...
     '\sigma^2=',th1(1),' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),...
     ' m; ',params.NyNx(1).*params.dydx(1),' m x ',...
     params.NyNx(2).*params.dydx(2),'m; dy=',params.dydx(1),'m, dx=',...
     params.dydx(2),'m; ',upper(masknm)),... 
    'Interpreter','latex');
   movev(ah,-0.03)                                                              
                                                                               
   % Pause to confirm the figure is nice before saving
   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo8_missing_Ny%iNx%i_dy%idx%i_s2%inu%irh%i_%s_1.eps',...
       params.NyNx,params.dydx,th1.*[1 10 1],date),'epsc')

   % Let's take a quick look at each taper 
   figure()
   colormap gray
   numrws=ceil(sqrt(mostmiss));
   [ah,ha,H]=krijetem(subnum(numrws,ceil(mostmiss./numrws)));
   for ind=1:size(ha,2)
     if ind<=mostmiss
       axes(ha(ind))
       imagesc(tapersv(:,:,ind))
       caxis([0 1])
       title(sprintf('%0.1f percent missing',100*prcmiss(ind)))
       longticks
       axis image
       shrink(ha(ind),1.02,1.02)
     else
       delete(ha(ind))
     end
   end
   movev([ha(1:mostmiss)],-0.03)
   sti(2)=sgtitle(...
     sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s%i%s%i%s%s',...
     '\sigma^2=',th1(1),' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),...
     ' m; ',params.NyNx(1).*params.dydx(1),' m x ',...
     params.NyNx(2).*params.dydx(2),'m; dy=',params.dydx(1),'m, dx=',...
     params.dydx(2),'m; ',upper(masknm)),...
    'Interpreter','latex');
   % Pause to confirm the figure is nice before saving
   keyboard
   saveas(gcf,...                                                               
     sprintf('covgammiosl_demo8_missing_Ny%iNx%i_dy%idx%i_s2%inu%irh%i_%s_2.eps',...
       params.NyNx,params.dydx,th1.*[1 10 1],date),'epsc') 
end
