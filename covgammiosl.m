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
% params.taper=ones(params.NyNx);
% params.taper(randi(prod(params.NyNx),floor(0.3*prod(params.NyNx)),1))=0;
% covgammiosl('demo4',[],params)
% % [IN PROGRESS] Demonstration of indexing for the DIAGONALS method:
% covgammiosl('demo5')
% % Study of the dependence of parameter covariance on grid-size (NyNx):
% covgammiosl('demo6')
% % Study of the dependence of parameter covariance on density (dydx):
% covgammiosl('demo7')
% % Study of the dependence of parameter covariance on missingness (taper),
% % either as random deletions or a growing masked region:
% covgammiosl('demo8')
% % Growing domain study comparing covthosl with the fisher inverse 
% covgammiosl('demo9')
%
% Last modified by fjsimons-at-princeton.edu, 04/08/2025
% Last modified by olwalbert-at-princeton.edu, 09/07/2026

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
  if isinf(params.blurs);params.blurs=-1;end
  Sbar=maternosp(th,params);

  % We also need blurred partial derivatives of the expected periodogram, which
  % we form from BLUROSY through MATERNOSY; initialize as [prod(NyNx) np] for
  % storage of the partials as vectors for each of the np parameters
  dSbardth=zeros([size(Sbar) np]);
  % One could use MAOSL but that is just using BLUROSY anyway
  for thnd=1:np
    if ifinv(thnd)
        % Only calculate partials for the parameters requested by ifinv
        dSbardthnd=blurosy(th,params,xver,'ef',[],thnd);
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
      numreals=500; 
      % Pre-allocate space for each score
      kgrads=zeros(numreals,np);

      % Simulate repetitively for the same parameters using circulant embedding,
      % returning the empirical periodogram of each realization, Semp
      runinpar=1;
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

      if isfield(params,'kiso') & ~isnan(params.kiso)
        % In the case that we have truncated the wavenumbers by setting
        % params.kiso, we will need to adjust the scaling. I should check that
        % this is accounted for in blurosy for Sbar and so on. If it is, this
        % becomes a counting issue only and what is below should be fine
        diferm(numel(kpp(kpp<=params.kiso))/numel(kpp)-kselbias(kpp,params))
        %lk=length(kpp(~~kpp)).*kselbias(kpp,params);
      % I think that the following commented out block will be necessary, but
      % not sufficient unless we build this same constraint into fishiosl which
      % is called by varianceofestimates, aka covthosl
      %  kgrads=deal(squeeze(sum((Sbar(kpp<=params.kiso&~~kpp)-...
      %                           Semps(kpp<=params.kiso&~~kpp,:)).*...
      %                        dSbardth(kpp<=params.kiso&~~kpp,:,:)./...
      %                        Sbar(kpp<=params.kiso&~~kpp).^2))); 
        lk=length(kpp(kpp<=params.kiso&~~kpp));
      else
        lk=length(kpp(~~kpp));
      end

      % The fair comparison would be between, essentially identical
      % diferm(ggrads,kgrads/lk)
      % [cov(ggrads) ; cov(kgrads/lk]
      % But in GAMMIOSL the zero wavenumber was excluded but in here not
      
      % Normalization factors in the case of masked and unit tapers
      normfact=lk^2;

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

      % The second term of the covariance of the periodogram via Isserlis' rule 
      % (eq. 9 of VoWE; also see Walden+1994)
      % In 1-D, EJJt=Um'*Cm*Um.'./NyNx;
      % In 2-D, EJJt=(Um*(Um*Cmn*Un.')*Un)'/prod(NyNx)
      % FIRST term in right hand side of A57 "transpose"
      EJJt_out=tensorprod(Um,tensorprod(Un,EJJht_in,1,2),1,2);
      EJJt_out=reshape(permute(EJJt_out,[3 1 2]),[prod(NyNx) prod(NyNx)]);
      
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
      normfact=prod(NyNx).^2;

      if isfield(params,'kiso')
        % [NOTE COPIED FROM CASE 1, IMPLEMENTATION NOT QUITE] 
        % In the case that we have truncated the wavenumbers by setting
        % params.kiso, we will need to adjust the scaling. I should check that
        % this is accounted for in blurosy for Sbar and so on. If it is, this
        % becomes a counting issue only and what is below should be fine
        kor=knums(params);
        normfact=normfact.*kselbias(kor,params).^2;
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
                                  max([1,mx+1]):min([Nx+mx,Nx]));
                  cdd_l2=abs(cdd{rnd}).^2;                                                  
                  cdd_l2=cdd_l2(:);                                                       
                  
                  % Use the same framework to select the antidiagonals from BLUROSY
                  fftcov_anti=blurosy(th,params,xver,'efd',[my2 mx2]);
                  fftcov_anti=v2s(fftcov_anti,params)*nf;
                  cad{rnd}=fftcov_anti(max([1,my2-Ny+2]):min([my2+1,Ny]),...
                                  max([1,mx2-Nx+2]):min([mx2+1,Nx]));
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
     normfact=prod(params.NyNx).^2;
     if isfield(params,'kiso')
       % [NOTE COPIED FROM CASE 1, IMPLEMENTATION NOT QUITE] 
       % In the case that we have truncated the wavenumbers by setting
       % params.kiso, we will need to adjust the scaling. I should check that
       % this is accounted for in blurosy for Sbar and so on. If it is, this
       % becomes a counting issue only and what is below should be fine
       kor=knums(params);
       normfact=normfact.*kselbias(kor,params).^2;
     end

     covg=s./normfact;
     % wkspc=whos; sum([wkspc.bytes])
  case 4
    % Attempt to clean up case 3 with simplifications 
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
    [x,y]=meshgrid(0:Nx-1,0:Ny-1);
    s_1=0;
    s_2=0;
    rnd=0;
    for my=(-Ny+1:Ny-1)
      dispprog=0;
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

          % Not the same [my my2 mx mx2]
          % This indexing is the same in test cases, but we should have it let
          % us know if this is ever not the case
         % diferm(1-prod([max([1,my+1])==max([1,my2-Ny+2])...
         %                min([Ny+my,Ny])==min([my2+1,Ny])...
         %                max([1,mx+1])==max([1,mx2-Nx+2])...
         %                min([Nx+mx,Nx])==min([mx2+1,Nx])]))

          subindy=[max([1,my+1]):min([Ny+my,Ny])];
          subindx=[max([1,mx+1]):min([Nx+mx,Nx])];
          % We acquire the offset FFT of the autocovariance via BLUROSY, which
          % we now provide method 'efd' and indices to calculate the offsets
          % in frequency as the tsto input argument
          fftcov_diag    =v2s(blurosy(th,params,xver,'efd',[my mx]),params)*nf;
          fftcov_diagsub =fftcov_diag(subindy,subindx);
          fftcov_diagsper=abs(fftcov_diagsub(:)).^2;

          % Use the same framework to select the antidiagonals from BLUROSY
          fftcov_anti    =v2s(blurosy(th,params,xver,'efd',[my2 mx2]),params)*nf;
          fftcov_antisub =fftcov_anti(subindy,subindx);
          fftcov_antisper=abs(fftcov_antisub(:)).^2;

          if normalize
              seq_ep =Sbar_n(ep_ind1).*Sbar_n(ep_ind2);
              ep_n3=ep_n(ep_ind3);
              seq_ep2=ep_n3.*flip(ep_n3,1);
          end

          % The partials indexed for the current `diagonal'
          f1=fax(ep_ind1,:); 
          f2=fax(ep_ind2,:);
          f3=flipud(f2);
          % Product of all partials in row-order:
          % [s2s2 s2nu s2rh nus2 nunu nurh rhs2 rhnu rhrh]
          f_seq=[f1(:,1).*f2 f1(:,2).*f2 f1(:,3).*f2];
          test=f1*f2';test(:);difer(f_seq-test)
          s_1=s_1+sum(fftcov_diagsper.*f_seq./seq_ep,1);

          f_seq2=[f2(:,1).*f3 f2(:,2).*f3 f2(:,3).*f3];
          s_2=s_2+sum(fftcov_antisper.*f_seq2./seq_ep2,1);
      end
    end
    s=reshape(s_1+s_2,np,np);
    normfact=prod(params.NyNx).^2;
    covg=s./normfact;
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
   % If the grid size is too large, we cannot use the dftmtx method. Set limit
   glim=120^2;

   try
     % Load a previous simulation suite
     datum=date;%'18-Feb-2025';
     % mleosl('demo2',datum)
     [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(datum,100);
     if isinf(p.blurs); p.blurs=-1; end
     if ~isfield(p,'taper'); p.taper=1; end
   catch
     try
       % Load a mat file of just thhats from working directory 
       % p.NyNx=[117 124];p.dydx=[1e3 1e3];p.blurs=Inf;th0=[1e2 1.5 2.8e3];
       % p.taper=1;
       % % Calculations for a full observed field:
       % pth='FigsFull/'
       % p.mask='full'; datum='30-Apr-2025';
       % % Calculations for the structured field with 66.7% observed:
       % pth='FigsStructure_66_7prc/'
       % p.mask='structure'; datum='10-May-2025';% '30-Apr-2025';
       % Calculations for the random field with 66.7% observed:
       % pth='FigsRandomDltn_66_7prc/'
       % p.mask='random'; datum='09-May-2025'; % '30-Apr-2025';
       % % Calculations for France with 66.7% observed:
       % pth='FigsFrance_66_7prc/'
       % p.mask='france'; datum='09-May-2025';
       % % Calculations for doves with 66.7% observed:
       % pth='FigsDoves_66_7prc/';
       % p.mask='doves'; datum='12-May-2025';
       % % Calculations for New Jersey with 66.7% observed:
       pth='FigsNewJersey_66_7prc/';
       p.mask='newjersey'; datum='12-May-2025';
       % Load everything we will need for regenerating the figures
       load(sprintf('%scovgammiosl_demo1vars_%s_%s.mat',pth,p.mask,datum))
     catch
       % Or calculate a new one
       % To match the paper example with blurosy demo 2,
       p.NyNx=[117 124];p.dydx=[1e3 1e3];p.blurs=Inf;
       th0=[1e2 1.5 2.8e3];

       if ~isfield(p,'mask')
         p.taper=1;
       elseif strcmp(p.mask,'france')
         scl=1.197; % corresponds to 66.7% data observed for bounded france on this grid size
         [~,I]=maskit(rand(p.NyNx),p,scl);
         p.taper=I;
       elseif strcmp(p.mask,'random')
         scl=0.667; % 66.7% data observed
         [~,I]=muckit(rand(p.NyNx),p,scl);
         p.taper=I;
       elseif strcmp(p.mask,'structure')
         I=load('MYMASKS/gebcostruct.mat'); % 66.7% data observed
         I=double(I.gebcostruct);
         p.taper=I;
       elseif strcmp(p.mask,'newjersey')
         scl=1.445; % 66.7% data observed
         [~,I]=maskit(rand(p.NyNx),p,scl);
         p.taper=I;
       elseif strcmp(p.mask,'doves')
         scl=2.03; % 66.7% data observed
         [~,I]=maskit(rand(p.NyNx),p,scl);
         p.taper=flipud(fliplr(I));
       end
       numreals=200;
       simp=p;simp.blurs=Inf;
       p.blurs=-1;
       for mnd=1:numreals
         thhat=NaN;
         while isnan(thhat)
           Hx=simulosl(th0,simp,0);
           [thhat,~,~,scl]=mleosl(Hx,[],p,[],[],th0,[1 1 1],0);
         end
         thhats(mnd,:)=thhat.*scl;
       end
       % Calculate parameter covariance for each method, storing the cross-parameter
       % covariance as unique variables
       covg1=covgammiosl(th0,p,1,ifinv);
       cov1=varianceofestimates(th0,p,covg1,ifinv);

       if prod(p.NyNx)<glim
         covg2=covgammiosl(th0,p,2,ifinv);
         cov2=varianceofestimates(th0,p,covg2,ifinv);
         cov2s2nu=matslice(cov2,[1 1 0]);
         cov2s2rh=matslice(cov2,[1 0 1]);
         cov2nurh=matslice(cov2,[0 1 1]);
       else
         cov2=NaN;
       end

       covg3=covgammiosl(th0,p,3,ifinv);
       cov3=varianceofestimates(th0,p,covg3,ifinv);

       % Unbiased outlier removal
       mdout=0;
       if mdout
         thhats=madout(thhats,1);
       end
       cove=cov(thhats);

       % Save
       if isfield(p,'mask')
         save(sprintf('%scovgammiosl_demo1vars_%s_%s.mat',pth,p.mask,date),...
              'thhats','p','cove','cov1','cov2','cov3','th0')
       else
         save(sprintf('covgammiosl_demo1vars_full_%s.mat',date),...
              'thhats','p','cove','cov1','cov2','cov3','th0')
       end
     end
   end
   % Slice it up
   cov1s2nu=matslice(cov1,[1 1 0]);
   cov1s2rh=matslice(cov1,[1 0 1]);
   cov1nurh=matslice(cov1,[0 1 1]);
   cov3s2nu=matslice(cov3,[1 1 0]);
   cov3s2rh=matslice(cov3,[1 0 1]);
   cov3nurh=matslice(cov3,[0 1 1]);
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
   
   if prod(p.NyNx)<glim
     % If the grid size is too large, we cannot use the dftmtx method
     [V,D]=eig(cov2s2nu);
     a21=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
     [V,D]=eig(cov2s2rh);
     a22=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
     [V,D]=eig(cov2nurh);
     a23=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
   end
   
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
   % for consistency between experiments that we are studying, set same limits
   shxlim={[20 200],[20 200],[0.8 2.6]};
   shylim={[0.8 2.6],[1900 3800],[1900 3800]};
   % Plot of s2 vs nu
   subplot(131)
   plot(thhats(:,1),thhats(:,2),'ko')
   hold on
   e1(1)=plot(a11(1,:)+mobs1(1),a11(2,:)+mobs1(2),'Color',clrs(1),'LineWidth',1);
   if prod(p.NyNx)<glim
     % If the grid size is too large, we cannot use the dftmtx method
     e2(1)=plot(a21(1,:)+mobs1(1),a21(2,:)+mobs1(2),'Color',clrs(2),'LineWidth',1);
   end
   e3(1)=plot(a31(1,:)+mobs1(1),a31(2,:)+mobs1(2),'Color',clrs(3),'LineWidth',1);
   ee(1)=plot(ae1(1,:)+mobs1(1),ae1(2,:)+mobs1(2),'Color','k','LineWidth',1);
   longticks
   xlabel('variance, $\sigma^2$','Interpreter','latex')
   ylabel('smoothness, $\nu$','Interpreter','latex')
   xlim(shxlim{1})
   ylim(shylim{1})

   % Plot of s2 vs rho 
   subplot(132)
   plot(thhats(:,1),thhats(:,3),'ko')
   hold on
   e1(2)=plot(a12(1,:)+mobs2(1),a12(2,:)+mobs2(2),'Color',clrs(1),'LineWidth',1);
   if prod(p.NyNx)<glim
     % If the grid size is too large, we cannot use the dftmtx method
     e2(2)=plot(a22(1,:)+mobs2(1),a22(2,:)+mobs2(2),'Color',clrs(2),'LineWidth',1);
   end
   e3(2)=plot(a32(1,:)+mobs2(1),a32(2,:)+mobs2(2),'Color',clrs(3),'LineWidth',1);
   ee(2)=plot(ae2(1,:)+mobs2(1),ae2(2,:)+mobs2(2),'Color','k','LineWidth',1);
   longticks
   xlabel('variance, $\sigma^2$','Interpreter','latex')
   ylabel('range, $\rho$','Interpreter','latex')
   xlim(shxlim{2})
   ylim(shylim{2})

   % Plot of nu vs rho 
   subplot(133)
   plot(thhats(:,2),thhats(:,3),'ko')
   hold on
   e1(3)=plot(a13(1,:)+mobs3(1),a13(2,:)+mobs3(2),'Color',clrs(1),'LineWidth',1);
   if prod(p.NyNx)<glim
     % If the grid size is too large, we cannot use the dftmtx method
     e2(3)=plot(a23(1,:)+mobs3(1),a23(2,:)+mobs3(2),'Color',clrs(2),'LineWidth',1);
   end
   e3(3)=plot(a33(1,:)+mobs3(1),a33(2,:)+mobs3(2),'Color',clrs(3),'LineWidth',1);
   ee(3)=plot(ae3(1,:)+mobs3(1),ae3(2,:)+mobs3(2),'Color','k','LineWidth',1);
   longticks
   xlabel('smoothness, $\nu$','Interpreter','latex')
   ylabel('range, $\rho$','Interpreter','latex')
   xlim(shxlim{3})
   ylim(shylim{3})

   if prod(p.NyNx)<glim
     % If the grid size is too large, we cannot use the dftmtx method
     [leg1,legic]=legend(sprintf('%i \\theta',size(thhats,1)),'sample','dftmtx',...
            'diagonals','empirical','box','off');
   else
     % Only 3 methods to report
     [leg1,legic]=legend(sprintf('%i \\theta',size(thhats,1)),'sample',...
            'diagonals','empirical','box','off');
   end
   % Find the 'line' objects
   movelcon(legic)

   % Resize the marker in the legend
   %set(legicls,'MarkerSize',20);
   sti=sgtitle(sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s',...
      '\sigma^2=',th0(1),' \mathrm{[unit]}^2, \nu=',th0(2),', \rho=',th0(3),...
      ' m; ',p.NyNx(1).*p.dydx(1),' m x ',p.NyNx(2).*p.dydx(2),'m'),... 
      'Interpreter','latex');
   % Save figure
   if isvarname('pth')
     saveas(gcf,...
       sprintf('%scovgammiosl_demo1_3pcrossplots_%ix%i_s2%inu%irh%i_%s_%s.eps',...
         pth,p.NyNx.*p.dydx,th0.*[1 10 1],p.mask,date),'epsc')
   else 
     saveas(gcf,...
       sprintf('covgammiosl_demo1_3pcrossplots_%ix%i_s2%inu%irh%i_%s_%s.eps',...
         p.NyNx.*p.dydx,th0.*[1 10 1],p.mask,date),'epsc')
   end
   
   % Plot the normalized covariance heatmaps for number of cov calc methods
   ncovhm=2;
keyboard
   figure()
   kelicol
   cmrange=[-1 1];
   labs={'s2','nu','rho'};
   labs=matslice(labs,ifinv);
   np=sum(ifinv);
   if ncovhm>2
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
     im(3)=imagesc(cov3./[diag(sqrt(cov3))*diag(sqrt(cov3))'],cmrange);
     title('diagonals')
     axis image

     % The dftmtx covariance
     if prod(p.NyNx)<glim
       % If the grid size is too large, we cannot use the dftmtx method
       ah(4)=subplot(224);
       im(4)=imagesc(cov2./[diag(sqrt(cov2))*diag(sqrt(cov2))'],cmrange);
       title('dftmtx')
       axis image
     end
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
     if prod(p.NyNx)<glim
       ah(4).Position(1)=(ah(2).Position(1)+ah(4).Position(1))/2+0.001;
     else
       ah(3).Position(1)=(ah(1).Position(1)+ah(3).Position(1))/2+0.001;
     end
   else
     % The empirical covariance
     ah(1)=subplot(121);
     im(1)=imagesc(matslice(cove,ifinv)./...
           [diag(sqrt(matslice(cove,ifinv)))*...
            diag(sqrt(matslice(cove,ifinv)))'],...
            cmrange);
     title('empirical')
     axis image
     % The diagonals covariance
     ah(2)=subplot(122);
     im(2)=imagesc(cov3./[diag(sqrt(cov3))*diag(sqrt(cov3))'],cmrange);
     title('predicted')
     axis image

     % Labels, colorbar, positioning, title
     set(ah(:),'xtick',1:np,'XTickLabel',labs)
     set(ah(:),'ytick',1:np,'YTickLabel',labs)
     longticks(ah(:))
     cb=colorbar;
     cb.Label.String='normalized';
     ah(2).Position(3)=ah(1).Position(3);
     moveh([ah(:)],-0.05)
   end

   ti=sgtitle(sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s',...
      '\sigma^2=',th0(1),' \mathrm{[unit]}^2, \nu=',th0(2),', \rho=',th0(3),...
      ' m; ',p.NyNx(1).*p.dydx(1),' m x ',p.NyNx(2).*p.dydx(2),'m'),... 
      'Interpreter','latex');

   % Save figure
   if isvarname('pth')
     saveas(gcf,...
       sprintf('%scovgammiosl_demo1_3pcorrplots_%ix%i_s2%inu%irh%i_%s_%s.eps',...
         pth,p.NyNx.*p.dydx,th0.*[1 10 1],p.mask,date),'epsc')
   else 
     saveas(gcf,...
       sprintf('covgammiosl_demo1_3pcorrplots_%ix%i_s2%inu%irh%i_%s_%s.eps',...
         p.NyNx.*p.dydx,th0.*[1 10 1],p.mask,date),'epsc')
   end
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
    params=[];params.NyNx=[94 97];params.dydx=[1 1];
    params.blurs=-1;params.taper=1;
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
    simparams=params;simparams.blurs=Inf;
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
    cov2=varianceofestimates(th1,params,covg2,ifinv);
    cov3=varianceofestimates(th1,params,covg3,ifinv);

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
    cov2=varianceofestimates(th1,params,covg2,ifinv);
    cov2=trilos(cov2);
    cov3=varianceofestimates(th1,params,covg3,ifinv);
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
    % We have the option to provide TH in the PARAMS input slot and PARAMS in
    % the METHOD input slot; each will be appropriately renamed
    defval('params',[1 0.5 2]) 
    defval('method',[]) 
    defval('ifinv',[1 1 1])
    th=params;
    params=method;
    rnddltn=1;
    rnddltn=0;

    if isempty(params)
      params.NyNx=[9 9];
      params.dydx=[1 1];
      params.blurs=-1;
      if rnddltn
        params.taper=ones(prod(params.NyNx),1);
        % the total number of deletions we want
        dff =floor(prod(params.NyNx)*1/3);
        % how many deletions we are lacking
        dffu=prod(params.NyNx)-sum(sum(params.taper))-dff;
        % randomly delete until we make a sufficient number of deletions
        % randi does not provide unique draws, so we loop
        while dffu<0
          params.taper(...
            randi(prod(params.NyNx),-dffu,1))=deal(0);
          dffu=prod(params.NyNx)-sum(sum(params.taper))-dff;
        end
        % reshape
        params.taper=v2s(params.taper,params);
      else
        params.taper=ones(params.NyNx);
      end
    end
  
    dydx=params.dydx;NyNx=params.NyNx;
    ys=0:NyNx(1)-1;xs=0:NyNx(2)-1;
    [grdx,grdy]=meshgrid(xs(:),ys(:));
    grdx=grdx(:).*dydx(2);
    grdy=grdy(:).*dydx(1);
    lagx=grdx-grdx';
    lagy=grdy-grdy';
    dxxp=sqrt(lagx.^2+lagy.^2);
    Cmn=maternosy(dxxp,th);
    if numel(params.taper)>1
      Cmn=Cmn.*params.taper(:).*params.taper(:)';
    end
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
    % These are not the same...
    diferm(EJJht_out-EJJht_out.') 
    % ...but these are
    diferm(EJJht_out-EJJht_out') 
    % This will be going to the "pSbar"
    EJJt_out=reshape(permute(...
               tensorprod(Um,tensorprod(Un,EJJht_in,1,2),1,2),...
             [3 1 2]),[prod(NyNx) prod(NyNx)]);
    Hk2cov=(abs(EJJht_out).^2+abs(EJJt_out).^2)./sum(sum(params.taper)).^2;
  
    % We will have two versions of this plot: the first (pltver==1) preserves
    % the original display of the terms involved in the Isserlis calculation,
    % the second (pltver==2) is for publication
    pltver=2;

    % While the first and second terms below sum to similar quantities, we can
    % quickly see wave-number dependent distinctions:
    figure(1)
    clf 
    figip=get(gcf,'InnerPosition');
    if pltver==1
      set(gcf,'InnerPosition',[figip(1:2) 1145 546])
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

      % Isserlis' periodogram covariance as the sum of the previous two terms 
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
    elseif pltver==2
      set(gcf,'InnerPosition',[figip(1:2) 1145 620])
      fs=11;
      ah=krijetem(subnum(2,4));
      axes(ah(1))
      imagesc(Cmn); caxis([0 1]); 
      cb(1)=colorbar;  cb(1).Label.String='[units^2]'; 
      cb(1).Location='south';
      cb(1).Ticks=0:0.2:1;
      cb(1).Position(2)=0.48;
      cb(1).Position(4)=0.015;
      cb(1).Position(3)=0.1436;
      cb(1).AxisLocation='out';
      cb(1).TickDirection='out';

      t(1)=title(append(repmat('\ ',1,17),...
           '$\mathcal{C}_{\boldmath{\theta}}(||\mathbf{x}-\mathbf{x''}||)$',...
           repmat('\ ',1,13),'$\rightarrow$'),...
           'FontSize',fs);
      t(1).Interpreter='latex';
      ah(1).XTick=[1-0.5 prod(params.NyNx)+0.5];
      ah(1).XTickLabel=[1 prod(params.NyNx)];
      ah(1).YTick=[1-0.5 prod(params.NyNx)+0.5];
      ah(1).YTickLabel=[1 prod(params.NyNx)];
      xlabel('$\mathbf{x''(:)}$','Interpreter','latex','FontSize',fs);
      ylabel('$\mathbf{x(:)}$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      axes(ah(2))
      imagesc(log(Hk2cov));
      cax=caxis; caxis(cax(2)-[10 0])
      cb(2)=colorbar; cb(2).Label.String='log, [units^4 length^2]';
      cb(2).Location='south';
      cb(2).Ticks=linspace(min(caxis),max(caxis),5);
      cb(2).TickLabels=round(cb(2).Ticks,1);
      cb(2).Position(2)=0.48;
      cb(2).Position(4)=0.015;
      cb(2).Position(3)=0.1436;
      cb(2).AxisLocation='out';
      cb(2).TickDirection='out';
      t(2)=title(append(repmat('\ ',1,10),...
                 '$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k''})|^2\}$',...
                 repmat('\ ',1,9),'$=$'),...
                 'FontSize',fs);
      t(2).Interpreter='latex';
      ah(2).XTick=[1 prod(params.NyNx)];
      ah(2).YTick=[1 prod(params.NyNx)];
      xlabel('$\mathbf{k''(:)}$','Interpreter','latex','FontSize',fs);
      ylabel('$\mathbf{k(:)}$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      axes(ah(3))
      imagesc(log(abs(EJJt_out./prod(NyNx)).^2));
      cax=caxis; caxis(cax(2)-[10 0])
      cb(3)=colorbar; cb(3).Label.String='log, [units^4 length^2]';
      cb(3).Location='south';
      cb(3).Ticks=linspace(min(caxis),max(caxis),5);
      cb(3).TickLabels=round(cb(3).Ticks,1);
      cb(3).Position(2)=0.48;
      cb(3).Position(4)=0.015;
      cb(3).Position(3)=0.1436;
      cb(3).AxisLocation='out';
      cb(3).TickDirection='out';
      t(3)=title(append(repmat('\ ',1,10),...
                 '$|\mathrm{cov}\{H(\mathbf{k}),H^*(\mathbf{k''})\}|^2$',...
                 repmat('\ ',1,10),'$+$'),...
                 'FontSize',fs);
      t(3).Interpreter='latex';
      ah(3).XTick=[1 prod(params.NyNx)];
      ah(3).YTick=[1 prod(params.NyNx)];
      xlabel('$\mathbf{k''(:)}$','Interpreter','latex','FontSize',fs);
      ylabel('$\mathbf{k(:)}$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      axes(ah(4))
      imagesc(log(abs(EJJht_out./prod(NyNx)).^2));
      cax=caxis; caxis(cax(2)-[10 0])
      cb(4)=colorbar; cb(4).Label.String='log, [units^4 length^2]';
      cb(4).Location='south';
      cb(4).Ticks=linspace(min(caxis),max(caxis),5);
      cb(4).TickLabels=round(cb(4).Ticks,1);
      cb(4).Position(2)=0.48;
      cb(4).Position(4)=0.015;
      cb(4).Position(3)=0.1436;
      cb(4).AxisLocation='out';
      cb(4).TickDirection='out';
      t(4)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k''})\}|^2$','FontSize',fs);
      t(4).Interpreter='latex';
      ah(4).XTick=[1 prod(params.NyNx)];
      ah(4).YTick=[1 prod(params.NyNx)];
      xlabel('$\mathbf{k''(:)}$','Interpreter','latex','FontSize',fs);
      ylabel('$\mathbf{k(:)}$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image
    end

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
    difer(diag(EJJht_out)./sum(sum(params.taper))-Sbar)

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
    difer(diag(EJJt_out)./sum(sum(params.taper))-pSbar)

    % The diagonal of the sum of the two predicted terms is the diagonal of A
    pSbar2=abs(pSbar).^2;

    diagpgmsum=v2s(Sbar2+pSbar2);
    difer(diag(A)-diagpgmsum(:))
    % Bottom line: Hk2cov is the sum of a covariance of Fourier coefficients and
    % a pseudocovariance of Fourier coefficients. The variances, aka their
    % diagonals, can be predicted from the blurred spectral density.

    if pltver==1
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
      xlabel('$k_x$','Interpreter','latex'); 
      ylabel('$k_y$','Interpreter','latex'); 
      longticks
      axis image

      ah(2)=subplot(222); imagesc(decibel(fftshift(v2s(Sbar2))))
      caxis([-40 0])
      t(2)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k})\}|^2$');
      t(2).Interpreter='latex';
      % Ticks represent indices only
      ah(2).XTick=[1 params.NyNx(1)];
      ah(2).YTick=[1 params.NyNx(2)];
      xlabel('$k_x$','Interpreter','latex'); 
      ylabel('$k_y$','Interpreter','latex'); 
      longticks
      axis image

      a=fftshift(v2s(Sbar2+pSbar2));
      ah(3)=subplot(223); imagesc(decibel(a));
      caxis([-40 0])
      t(3)=title('$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k})|^2\}$');
      t(3).Interpreter='latex';
      ah(3).XTick=[1 params.NyNx(1)];
      ah(3).YTick=[1 params.NyNx(2)];
      xlabel('$k_x$','Interpreter','latex'); 
      ylabel('$k_y$','Interpreter','latex'); 
      longticks
      axis image

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
      xlabel('$k_x$','Interpreter','latex'); 
      ylabel('$k_y$','Interpreter','latex'); 
      longticks
      axis image

      % layout(ah(3),0.5,'m','x')
      
      % % Check if Hermitian
      % Sbar2=fftshift(v2s(Sbar2));
      % pSbar2=fftshift(v2s(pSbar2));
      % hermcheck(Sbar2)
      % hermcheck(pSbar2)
      % % Check if positive definite 
      % flg=defcheck(Sbar2); diferm(flg>1)
      % pflg=defcheck(pSbar2); diferm(pflg>1)

      saveas(gcf,...
        sprintf('covgammiosl_demo4_dftmtx_%ix%i_s2%inu%irh%i_%s.eps',...
          params.NyNx.*params.dydx,th.*[1 10 1],date),'epsc')
   elseif pltver==2
      axes(ah(5))
      b=v2s(blurosy(th,params)*(2*pi)^2).^2;
      [~,i]=max(b(:)); b(i)=b(i)*2;
      imagesc(decibel(b))
      caxis([-40 0])
      cb(5)=colorbar; cb(5).Label.String='log, [units^4 length^2]';
                        %repmat('\ ',1,21),...
                        %repmat('\ ',1,18),'$=$'),...
      cb(5).Location='south';
      cb(5).Ticks=linspace(min(caxis),max(caxis),5);
      cb(5).TickLabels=round(cb(5).Ticks,1);
      cb(5).Position(2)=0.01;
      cb(5).Position(4)=0.015;
      cb(5).Position(3)=0.1436;
      cb(5).AxisLocation='out';
      cb(5).TickDirection='out';
      t(5)=title('$\bar{\mathcal{S}_{\boldmath{\theta}}}^2(\mathbf{k})$',...
                 'FontSize',fs);
      t(5).Interpreter='latex';
      ah(5).XTick=[1-0.5 params.NyNx(1)+0.5];
      ah(5).XTickLabel=[1 params.NyNx(1)];
      ah(5).YTick=[1-0.5 params.NyNx(2)+0.5];
      ah(5).YTickLabel=[1 params.NyNx(2)];
      xlabel('$k_x$','Interpreter','latex','FontSize',fs); 
      ylabel('$k_y$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      axes(ah(6))
      a=fftshift(v2s(Sbar2+pSbar2));
      imagesc(decibel(a));
      caxis([-40 0])
      cb(6)=colorbar; cb(6).Label.String='log, [units^4 length^2]';
      cb(6).Location='south';
      cb(6).Ticks=linspace(min(caxis),max(caxis),5);
      cb(6).TickLabels=round(cb(6).Ticks,1);
      cb(6).Position(2)=0.01;
      cb(6).Position(4)=0.015;
      cb(6).Position(3)=0.1436;
      cb(6).AxisLocation='out';
      cb(6).TickDirection='out';
      t(6)=title(append(repmat('\ ',1,9),...
                        '$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k})|^2\}$',...
                        repmat('\ ',1,10),'$=$'),...
                        'FontSize',fs);
      t(6).Interpreter='latex';
      ah(6).XTick=[1-0.5 params.NyNx(1)+0.5];
      ah(6).XTickLabel=[1 params.NyNx(1)];
      ah(6).YTick=[1-0.5 params.NyNx(2)+0.5];
      ah(6).YTickLabel=[1 params.NyNx(2)];
      xlabel('$k_x$','Interpreter','latex','FontSize',fs); 
      ylabel('$k_y$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      axes(ah(7))
      imagesc(decibel(fftshift(v2s(pSbar2))))
      caxis([-40 0])
      cb(7)=colorbar; cb(7).Label.String='log, [units^4 length^2]';
      cb(7).Location='south';
      cb(7).Ticks=linspace(min(caxis),max(caxis),5);
      cb(7).TickLabels=round(cb(7).Ticks,1);
      cb(7).Position(2)=0.01;
      cb(7).Position(4)=0.015;
      cb(7).Position(3)=0.1436;
      cb(7).AxisLocation='out';
      cb(7).TickDirection='out';
      t(7)=title(append(repmat('\ ',1,10),...
                        '$|\mathrm{cov}\{H(\mathbf{k}),H^*(\mathbf{k})\}|^2$',...
                        repmat('\ ',1,10),'$+$'),...
                        'FontSize',fs);
      t(7).Interpreter='latex';
      ah(7).XTick=[1-0.5 params.NyNx(1)+0.5];
      ah(7).XTickLabel=[1 params.NyNx(1)];
      ah(7).YTick=[1-0.5 params.NyNx(2)+0.5];
      ah(7).YTickLabel=[1 params.NyNx(2)];
      xlabel('$k_x$','Interpreter','latex','FontSize',fs); 
      ylabel('$k_y$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      axes(ah(8))
      imagesc(decibel(fftshift(v2s(Sbar2))))
      caxis([-40 0])
      cb(8)=colorbar; cb(8).Label.String='log, [units^4 length^2]';
      cb(8).Location='south';
      cb(8).Ticks=linspace(min(caxis),max(caxis),5);
      cb(8).TickLabels=round(cb(8).Ticks,1);
      cb(8).Position(2)=0.01;
      cb(8).Position(4)=0.015;
      cb(8).Position(3)=0.1436;
      cb(8).AxisLocation='out';
      cb(8).TickDirection='out';
      t(8)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k})\}|^2$','FontSize',fs);
      t(8).Interpreter='latex';
      % Ticks represent indices only
      ah(8).XTick=[1-0.5 params.NyNx(1)+0.5];
      ah(8).XTickLabel=[1 params.NyNx(1)];
      ah(8).YTick=[1-0.5 params.NyNx(2)+0.5];
      ah(8).YTickLabel=[1 params.NyNx(2)];
      xlabel('$k_x$','Interpreter','latex','FontSize',fs); 
      ylabel('$k_y$','Interpreter','latex','FontSize',fs); 
      longticks
      axis image

      movev(cb,0.055)
      movev(ah,0.01)
     % movev([ah(1:4) cb(1:4)],0.06)
     % movev([cb(1:4)],0.03)
     % movev([ah(5:8) cb(5:8)],0.02)
     % movev([cb(5:8)],0.03)

     % movev([t(1:4)],-1.3)
     % movev([t(5:8)],-0.2)
     %[bh,txh]=label(ah,'ul',9,0);
     set(gcf,'Renderer','painters')
keyboard
      saveas(gcf,...
        sprintf('covgammiosl_demo4_dftmtx_%ix%i_s2%inu%irh%i_%s_v2.eps',...
          params.NyNx.*params.dydx,th.*[1 10 1],date),'epsc')

      pltavg=1;                                                                 
      if pltavg                                                                 
        varns={Cmn,Hk2cov,...                                                   
               abs(EJJt_out./prod(NyNx)).^2,abs(EJJht_out./prod(NyNx)).^2,...   
               b,fftshift(v2s(Sbar2+pSbar2)),...                                
               fftshift(v2s(pSbar2)),fftshift(v2s(Sbar2)),ah};                     
        varargout=varns(1:nargout);                                             
      end      
   end
elseif strcmp(th,'demo4b')                                                      
    N=250;                                                                      
    params=[];params.dydx=[1 1];params.blurs=-1;                                
    params.NyNx=[9 9];                                                          
    propdltn=1/3;                                                               
    pttrn='random';                                                             
    pttrn='newjersey';                                                         
    if strcmp(pttrn,'random')                                                   
      for ind=1:N                                                               
        tpr=ones(params.NyNx);                                                  
        tpr(randi(prod(params.NyNx),[floor(prod(params.NyNx)*propdltn) 1]))=deal(0);
        while sum(tpr(:))/81>1-propdltn                                         
          tpr(randi(prod(params.NyNx),1))=deal(0);                              
        end                                                                     
        params.taper=tpr;                                                       
        tprs{ind}=tpr;                                                          
        [a{ind},b{ind},c{ind},d{ind},e{ind},f{ind},g{ind},h{ind},ah{ind}]=...           
                covgammiosl('demo4',[],params);                                 
      end                                                                       
    elseif strcmp(pttrn,'newjersey')                                            
      params.mask=pttrn;
      [~,I]=maskit(ones(params.NyNx),params,[]);
      % params.taper=~I; %anti new jersey
      params.taper=I;
      ind=1;
      propdltn=sum(params.taper(:))/prod(params.NyNx);
      [a{ind},b{ind},c{ind},d{ind},e{ind},f{ind},g{ind},h{ind},ah{ind}]=...           
                covgammiosl('demo4',[],params);                                 
      keyboard
      [bh,txh]=label(ah{ind},'ul',9,4);
      set(gcf,'Renderer','painters')
    end                                                                         
    keyboard
    pers={a,b,c,d,e,f,g,h};                                                     
    for jnd=1:size(pers,2)                                                      
      avgpers{jnd}=mean(cat(3,pers{jnd}{:}),3);                                 
    end                                                                         
    % plot the averages using the same layout as pltver==2 from demo4           
    fs=11;                                                                      
    figure(2)                                                                   
    figip=get(gcf,'InnerPosition');
    set(gcf,'InnerPosition',[figip(1:2) 1145 546])
    ah=krijetem(subnum(2,4));                                                   
    axes(ah(1))                                                                 
    imagesc(avgpers{1}); caxis([0 1]);                                          
    t(1)=title(append(repmat('\ ',1,17),...                                     
         '$\mathcal{C}_{\boldmath{\theta}}(||\mathbf{x}-\mathbf{x''}||)$',...   
         repmat('\ ',1,13),'$\rightarrow$'),...                                 
         'FontSize',fs);                                                        
    t(1).Interpreter='latex';                                                   
    ah(1).XTick=[1-0.5 prod(params.NyNx)+0.5];                                  
    ah(1).XTickLabel=[1 prod(params.NyNx)];                                     
    ah(1).YTick=[1-0.5 prod(params.NyNx)+0.5];                                  
    ah(1).YTickLabel=[1 prod(params.NyNx)];                                     
    xlabel('$\mathbf{x''(:)}$','Interpreter','latex','FontSize',fs);            
    ylabel('$\mathbf{x(:)}$','Interpreter','latex','FontSize',fs);              
    longticks                                                                   
    axis image                                                                  
                                                                                
    axes(ah(2))                                                                 
    imagesc(log(avgpers{2}))                                                    
    cax=caxis; caxis(cax(2)-[10 0])                                             
    t(2)=title(append(repmat('\ ',1,10),...                                     
               '$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k''})|^2\}$',...    
               repmat('\ ',1,9),'$=$'),...                                      
               'FontSize',fs);                                                  
    t(2).Interpreter='latex';                                                   
    ah(2).XTick=[1 prod(params.NyNx)];                                          
    ah(2).YTick=[1 prod(params.NyNx)];                                          
    xlabel('$\mathbf{k''(:)}$','Interpreter','latex','FontSize',fs);            
    ylabel('$\mathbf{k(:)}$','Interpreter','latex','FontSize',fs);              
    longticks                                                                   
    axis image                                   

    axes(ah(3))
    imagesc(log(avgpers{3}))
    cax=caxis; caxis(cax(2)-[10 0])
    t(3)=title(append(repmat('\ ',1,10),...
               '$|\mathrm{cov}\{H(\mathbf{k}),H^*(\mathbf{k''})\}|^2$',...
               repmat('\ ',1,10),'$+$'),...
               'FontSize',fs);
    t(3).Interpreter='latex';
    ah(3).XTick=[1 prod(params.NyNx)];
    ah(3).YTick=[1 prod(params.NyNx)];
    xlabel('$\mathbf{k''(:)}$','Interpreter','latex','FontSize',fs);
    ylabel('$\mathbf{k(:)}$','Interpreter','latex','FontSize',fs);
    longticks
    axis image

    axes(ah(4))
    imagesc(log(avgpers{4}))
    cax=caxis; caxis(cax(2)-[10 0])
    t(4)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k''})\}|^2$','FontSize',fs);
    t(4).Interpreter='latex';
    ah(4).XTick=[1 prod(params.NyNx)];
    ah(4).YTick=[1 prod(params.NyNx)];
    xlabel('$\mathbf{k''(:)}$','Interpreter','latex','FontSize',fs);
    ylabel('$\mathbf{k(:)}$','Interpreter','latex','FontSize',fs);
    longticks
    axis image

    axes(ah(5))
    imagesc(decibel(avgpers{5}))
    caxis([-40 0])
    t(5)=title('$\bar{\mathcal{S}_{\boldmath{\theta}}}^2(\mathbf{k})$',...
               'FontSize',fs);
    t(5).Interpreter='latex';
    ah(5).XTick=[1-0.5 params.NyNx(1)+0.5];
    ah(5).XTickLabel=[1 params.NyNx(1)];
    ah(5).YTick=[1-0.5 params.NyNx(2)+0.5];
    ah(5).YTickLabel=[1 params.NyNx(2)];
    xlabel('$k_x$','Interpreter','latex','FontSize',fs);
    ylabel('$k_y$','Interpreter','latex','FontSize',fs);
    longticks
    axis image

    axes(ah(6))
    imagesc(decibel(avgpers{6}))
    caxis([-40 0])
    t(6)=title(append(repmat('\ ',1,9),...
                      '$\mathrm{cov}\{|H(\mathbf{k})|^2,|H(\mathbf{k})|^2\}$',...
                      repmat('\ ',1,10),'$=$'),...
                      'FontSize',fs);
    t(6).Interpreter='latex';
    ah(6).XTick=[1-0.5 params.NyNx(1)+0.5];
    ah(6).XTickLabel=[1 params.NyNx(1)];
    ah(6).YTick=[1-0.5 params.NyNx(2)+0.5];
    ah(6).YTickLabel=[1 params.NyNx(2)];
    xlabel('$k_x$','Interpreter','latex','FontSize',fs);
    ylabel('$k_y$','Interpreter','latex','FontSize',fs);
    longticks
    axis image

    axes(ah(7))
    imagesc(decibel(avgpers{7}))
    caxis([-40 0])
    t(7)=title(append(repmat('\ ',1,10),...
                      '$|\mathrm{cov}\{H(\mathbf{k}),H^*(\mathbf{k})\}|^2$',...
                      repmat('\ ',1,10),'$+$'),...
                      'FontSize',fs);
    t(7).Interpreter='latex';
    ah(7).XTick=[1-0.5 params.NyNx(1)+0.5];
    ah(7).XTickLabel=[1 params.NyNx(1)];
    ah(7).YTick=[1-0.5 params.NyNx(2)+0.5];
    ah(7).YTickLabel=[1 params.NyNx(2)];
    xlabel('$k_x$','Interpreter','latex','FontSize',fs);
    ylabel('$k_y$','Interpreter','latex','FontSize',fs);
    longticks
    axis image

    axes(ah(8))
    imagesc(decibel(avgpers{8}))
    caxis([-40 0])
    t(8)=title('$|\mathrm{cov}\{H(\mathbf{k}),H(\mathbf{k})\}|^2$','FontSize',fs);
    t(8).Interpreter='latex';
    % Ticks represent indices only
    ah(8).XTick=[1-0.5 params.NyNx(1)+0.5];
    ah(8).XTickLabel=[1 params.NyNx(1)];
    ah(8).YTick=[1-0.5 params.NyNx(2)+0.5];
    ah(8).YTickLabel=[1 params.NyNx(2)];
    xlabel('$k_x$','Interpreter','latex','FontSize',fs);
    ylabel('$k_y$','Interpreter','latex','FontSize',fs);
    longticks
    axis image
    movev([t(1:4)],-1.3)
    movev([t(5:8)],-0.2)

    saveas(gcf,...
        sprintf('covgammiosl_demo4b_%i_%ix%i_%s.eps',...
        round(100-propdltn*100),params.NyNx.*params.dydx,date),'epsc')
    keyboard

    mnd=7;[mx,el]=max(avgpers{mnd},[],2)
    figure(3)
    avgtprs=mean(cat(3,tprs{:}),3);
    imagesc(avgtprs{9})
    keyboard

elseif strcmp(th,'demo5')
   % [IN PROGRESS] Visual demonstration of the indexing scheme implemented in the DIAGONALS
   % method
   % TODO: COPY IN SCRATCH/OLIVIASJ.M
elseif strcmp(th,'demo6')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the GRID SIZE, which is presented in terms
   % of pi rho by pi rho. Plot the bias, (co)variance, and rmse from the 
   % empirical study and calculations, when applicable, for each parameter.
   params=[];params.dydx=[1 1];params.blurs=-1;params.taper=1;ifinv=[1 1 1];
   %numreals=200;
   numreals=100;
   np=3;

   % Labels
   labs1={'$\sigma^2$','$\nu$','$\rho$'};
   labs={labs1{:},'\sigma^2,\nu','\sigma^2,\rho', '\nu,\rho',labs1{:},labs1{:}};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];

   % For now, small values so that we have a chance at forming the data vector
   % and finding an estimate for small grids 
   rh=6;
   rh=4; % 01/20/26 -- this needs to be small enough for circulant embedding
   %rh=16;
   th1=[3 0.75 rh];
   % Sample growing grid length as multiples of pi rho
   glens=5:8;
   Ns=(2.^glens).^2;
   log2glens=sqrt(Ns);
   % Second set of vectors for the second dataset
   glens2=log2(sqrt(Ns*2/3));
   log2glens2=sqrt(Ns*2/3);

   thhats=zeros(numreals,3,size(glens,2));
   % for the full field: 
   %datum='30-Apr-2025';
   % for 'random' 66.7% of field: 
   %datum='11-May-2025'; params.mask='random';
   datum='21-Jan-2026'; params.mask='random';

   % Flag for paper figure that will show the field fully observed and with a
   % constant ratio of random deletions on the same plot
   flg=1;
   if flg==1
     % We will let the fully observed field be loaded as usual
     datum='30-Apr-2025';
     % We will load the field with random deletions here and denote its
     % variables with a '2'
     datum2='11-May-2025'; params2.mask='random';
     datum2='21-Jan-2026';params2.mask='random';
     demo6vars2=load(sprintf('covgammiosl_demo6_calcs_%s.mat',datum2));
     thhats2=demo6vars2.thhats;
    % cov1t2=demo6vars2.cov1t;
     cov3t2=demo6vars2.cov3t;
     covempt2=demo6vars2.covempt;
     Ns2=demo6vars2.Ns;
     % Distinguish each dataset by color
     clrs=["k","k","k"];
     clrs2=["#1B9E77","#1B9E77","#1B9E77"];
   end
   try
     % If we have already calculated the estimates, load them
    % for the full field: datum='20-Jan-2026'; needs more reals
    % for the randomly deleted 66.7% field: datum='21-Jan-2026';
     datum='20-Jan-2026';
     demo6vars=load(sprintf('covgammiosl_demo6_calcs_%s.mat',datum));
     thhats=demo6vars.thhats;
     %cov1t=demo6vars.cov1t;
     %try
     %  % We might not have been able to calculate this depending on the grid
     %  % size
     %  cov2t=demo6vars.cov2t;
     %  incdft=1;
     %catch
     %  incdft=0;
     %end
     cov3t=demo6vars.cov3t;
     covempt=demo6vars.covempt;
     Ns=demo6vars.Ns;
   catch
     for ind=1:size(log2glens,2)
       datum=date;
       datum='22-Jan-2026';params2.mask='random';
       % Create a grid whose length will increase logarithmically 
       params.NyNx=[log2glens(ind) log2glens(ind)];
       if isfield(params,'mask'); if strcmp(params.mask,'random'); 
         [~,I]=muckit(randn(params.NyNx),params,0.667); 
         tapersv{ind}=I;
         params.taper=tapersv{ind};
       end; end
       % Simulate via circulant embedding
       simparams=[];simparams.dydx=params.dydx;simparams.NyNx=params.NyNx;
       simparams.blurs=Inf;simparams.taper=1;
       for mnd=1:numreals
         thhat=NaN;
         while isnan(thhat)
           Hx=simulosl(th1,simparams,0);
           [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
         end
         thhats(mnd,:,ind)=thhat.*scl;
       end
     end
     mdout=0;
     if mdout
       % Apply MADOUT to robustly remove outliers; if outliers are removed,
       % continue simulating until acquiring NUMREALS estimates
       for ind=1:size(log2glens,2)
         % Create a grid whose length will increase logarithmically 
         params.NyNx=[log2glens(ind) log2glens(ind)];
         if isfield(params,'mask'); if strcmp(params.mask,'random'); 
             [~,I]=muckit(randn(params.NyNx),params,0.667); 
             tapersv{ind}=I;
             params.taper=tapersv{ind};
         end; end
         % Simulate via circulant embedding
         simparams=params;simparams.blurs=Inf;
         thhatst{ind}=madout(thhats(:,:,ind),1);
         while size(thhatst{ind},1)<numreals
           for jnd=1:(numreals-size(thhatst{ind},1))
             Hx=simulosl(th1,simparams,0);
             thhat=NaN;
             while isnan(thhat)
               [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],[],ifinv,0);
             end
             thhatspls(jnd,:)=thhat.*scl;
           end
           % Recreate thhats(:,:,ind) with the recalculated estimates
           thhats(:,:,ind)=[thhatst{ind};thhatspls];
           thhatspls=[];
           thhatst{ind}=madout(thhats(:,:,ind),1);
         end
       end
     end
     % Calculate parameter covariance up to 4 different ways. We may choose not to
     % include the dftmtx method in favor of looking at larger grids (the
     % diagonals method will provide identical results) 
     for ind=1:size(log2glens,2)
       % Create a grid whose length will increase logarithmically 
       params.NyNx=[log2glens(ind) log2glens(ind)];
       if isfield(params,'mask'); if strcmp(params.mask,'random'); 
           params.taper=tapersv{ind};
       end; end
       % Calculate the empirical covariance from the ensemble
       covemp=nancov(thhats(:,:,ind));
       % Calculate the approximat covariance with the 'sample' method
       %cov1c=covthosl(th1,params,1,ifinv);
       %if incdft
       %  % Calculate the exact covariance with the 'dftmtx' method
       %  cov2c=covthosl(th1,params,2,ifinv);
       %end
       % Calculate the exact covariance with the diagonals method
       meth=3;
       disp(sprintf('Calculating COVTHOSL for meth %i and NyNx %i %i',...
         meth,params.NyNx))
       cov3c=covthosl(th1,params,meth,ifinv);
       % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
       covempt(ind,:)=trilos(covemp);
       cov3t(ind,:)  =trilos(cov3c);
       % cov1t(ind,:)  =trilos(cov1c);
       % cov2t(ind,:)  =trilos(cov2c);
     end
     % Consider saving your calculations by setting 'svcalcs'
     svcalcs=1;
     if svcalcs
       save(sprintf('covgammiosl_demo6_calcs_%s.mat',datum),...
            'thhats','cov3t','covempt','params','th1','Ns')
     end
   end

   clf;
   figip=get(gcf,'InnerPosition');
   set(gcf,'InnerPosition',[figip(1:2) 950 950])
   [ah,ha,H]=krijetem(subnum(4,3));
   fs=11;
   xla='$\sqrt{K}$';
   gv=0.8*glens(1):0.1:1.2*glens(end);
   gv2=0.8*glens2(1):0.1:1.2*glens2(end);
   xtks=glens;
   xtklbs=string(log2glens);
   xxtklbs=string(ceil(log2glens/(pi*rh)));
   xxtks=log2(ceil(log2glens/(pi*rh))*pi*rh);
   % I used to be interested in showing that each of the covthosl calculations
   % are essentially identical, but since this is established, we won't bother
   % to muddy the plot:
   covana=cov3t;
   covana2=cov3t2;
   % Complete the calculations. Fit a slope to the covariance calculations as a
   % function of grid length
   for jnd=1:6
     mbfemp(jnd,:)     =polyfit(glens, 10*log10(abs(covempt(:,jnd))),1);
     mbfcovana(jnd,:)  =polyfit(glens, 10*log10(abs(covana(:,jnd))),1);
     predemp(jnd,:)    =polyval(mbfemp(jnd,:),   gv);
     predcovana(jnd,:) =polyval(mbfcovana(jnd,:),  gv);
     mbfemp2(jnd,:)    =polyfit(glens2,10*log10(abs(covempt2(:,jnd))),1);
     mbfcovana2(jnd,:) =polyfit(glens2,10*log10(abs(covana2(:,jnd))),1);
     predemp2(jnd,:)   =polyval(mbfemp2(jnd,:),   gv2);
     predcovana2(jnd,:)=polyval(mbfcovana2(jnd,:),gv2);
     if jnd<4
       % rmse
       rmseemp(jnd,:)   =sqrt(squeeze(var(thhats(:,jnd,:)))+...
                             squeeze(mean(thhats(:,jnd,:))-th1(jnd)).^2);
       mbfrmse(jnd,:)   =polyfit(glens,10*log10(abs(rmseemp(jnd,:))),1);
       predrmse(jnd,:)  =polyval(mbfrmse(jnd,:),   gv);
       rmseemp2(jnd,:)  =sqrt(squeeze(var(thhats2(:,jnd,:)))+...
                             squeeze(mean(thhats2(:,jnd,:))-th1(jnd)).^2);
       mbfrmse2(jnd,:)  =polyfit(glens2,10*log10(abs(rmseemp2(jnd,:))),1);
       predrmse2(jnd,:) =polyval(mbfrmse2(jnd,:),   gv2);
     end
   end

   for ind=1:12
     axes(ah(ind))
     hold on
     if ind<4
       % First row. Plot the relative bias of the estimates from the empirical study
       relbsme(:,ind)=squeeze(mean(thhats(:,ind,:)-th1(ind)))/th1(ind);
       relbssd(:,ind)=squeeze(sqrt(var(thhats(:,ind,:)-th1(ind),"omitnan"))/...
                      th1(ind));
       relbsme2(:,ind)=squeeze(mean(thhats2(:,ind,:)-th1(ind)))/th1(ind);
       relbssd2(:,ind)=squeeze(sqrt(var(thhats2(:,ind,:)-th1(ind),"omitnan"))/...
                      th1(ind));
       % First dataset (100% observed grid)
       % Errorbar in terms of the standard deviation
       pl1(ind)=errorbar(glens,relbsme(:,ind),relbssd(:,ind),'ko','CapSize',5,...
                   'DisplayName','100% obs');
       % Add the bias markers for the second dataset
       pl12(ind)=errorbar(glens2,relbsme2(:,ind),relbssd2(:,ind),'o',...
                   'MarkerEdgeColor',clrs2{1},'MarkerFaceColor',clrs2{1},...
                   'LineWidth',1,'CapSize',5,'Color',clrs2{1},...
                   'DisplayName','66.7% obs');
       xlim(minmax([glens glens2]).*[0.99 1.01])
      % We don't need this anymore; it affects the legend
      % if ind==2
      %   % the bias in nu for the full grid is much smaller, so lets pull it on
      %   % top
      %   uistack(pl1(2),'up',1)
      % end
       % Add the 0-line where the estimate equals the truth
       yline(0,'HandleVisibility','off')
       if ind==1
         ylabel('$(\langle\hat{\theta}\rangle-\theta_0)/\theta_0$','interpreter','latex');
      % elseif ind==2
      %   keyboard
      %   ylim([-0.2,max(nanmean(thhats(:,ind,:))-th1(ind))+0.2])
       end
       ti(ind)=title(labs{ind},'Interpreter','latex','Units','pixels');
       [leg1(ind),legic]=legend('interpreter','tex','EdgeColor','none',...
                            'NumColumns',2,'FontSize',7.5,...
                            'BackgroundAlpha',0,'box','off');

       % Have the legend markers take up less space
       movelcon(legic)

       if ind==3
         % Cosmetics for the first row, including second xaxis
         byl=max(max(abs(cell2mat(ylim(ah([1:3]))))));
         byl=[-byl byl];
         byl=[-0.8 0.8];
         for knd=1:ind
           axes(ah(knd))
           ylim([byl])
           longticks
           if knd>1
             set(gca,'YTickLabel','')
           end
           % Set the top x-axis
           xlb='Length [$\pi \rho$]';
           [axx(knd),xl(knd),yl(knd)]=...
             xtraxis(ah(knd),xxtks,xxtklbs,xlb);
           set(xl(knd),'Interpreter','latex')
           longticks
           ah(knd).XRuler.TickLabelGapOffset=-2;
           ah(knd).XTickLabelRotation=0;
           axx(knd).XMinorTick='off';
           axx(knd).XRuler.TickLabelGapOffset=-3;
         end
       end
     elseif ind<10
       % The first and second rows. Plot the absolute (co)variance.
       % The first dataset (100% observed grid)
       jnd=ind-np;
       pe{1,jnd}=plot(glens,10*log10(abs(covempt(:,jnd))),'ko',...
                       'DisplayName',sprintf('100%% obs %0.1f',abs(mbfemp(jnd,1))));
       pa{1,jnd}=plot(glens,10*log10(abs(covana(:,jnd))) ,'^','MarkerEdgeColor',clrs{1},...
              'MarkerFaceColor','None',...
              'DisplayName',sprintf('100%% pred %0.1f',abs(mbfcovana(jnd,1))));
       % Plot the absolute covariance calculated from the ensemble and covthosl
       ppe{1,jnd}=plot(gv,predemp(jnd,:),'Color','k','HandleVisibility','off');
       ppa{2,jnd}=plot(gv,predcovana(jnd,:),'Color',clrs{1},'HandleVisibility','off');

       % The second dataset (66.7% observed grid, random deletions)
       pe{2,jnd}=plot(glens2,10*log10(abs(covempt2(:,jnd))),'o','MarkerEdgeColor',clrs2{1},...
                      'MarkerFaceColor',clrs2{1},...
                      'DisplayName',sprintf('66.7%% obs %0.1f',abs(mbfemp2(jnd,1))));
       pa{2,jnd}=plot(glens2,10*log10(abs(covana2(:,jnd))),'^','MarkerEdgeColor',clrs2{1},...
                 'MarkerFaceColor',clrs2{1},...
                 'DisplayName',sprintf('66.7%% pred %0.1f',abs(mbfcovana2(jnd,1))));
       % Plot the absolute covariance calculated from the ensemble and covthosl
       ppe{1,jnd}=plot(gv2,predemp2(jnd,:),'Color',clrs2{3},'LineStyle','-',...
                       'HandleVisibility','off');
       ppa{2,jnd}=plot(gv2,predcovana2(jnd,:),'Color',clrs2{3},'LineStyle','-',...
                       'HandleVisibility','off');
       %set(gca,'YLim',[1e-5 1e7],'YTick',[1e-3 1e0 1e3 1e6])
       if ~ismember(jnd,[1 4])
         set(gca,'YTickLabel','')
       end
       [leg1(ind),legic]=legend('interpreter','tex','EdgeColor','none',...
                                'NumColumns',2,'FontSize',7.5,...
                                'BackgroundAlpha',0,'box','off');
       leg1(ind).AutoUpdate='off';
       % Have the legend markers take up less space
       movelcon(legic)
       
       % Labels, titles, ticks, common limits
       if jnd==1
          ylabel('$10\times \mathrm{log}_{10}(|\mathrm{var}\{\theta_i\}|)$',...
               'Interpreter','latex','FontSize',fs)
       elseif jnd==4
          ylabel('$10\times \mathrm{log}_{10}(|\mathrm{cov}\{\theta_i,\theta_j\}|$)',...
                 'Interpreter','latex','FontSize',fs)
       end
       if jnd==6
         % Set consistent axes
         allylims=[ylim(ah([np+1:np+6]))];
         allylims=minmax([allylims{:}]);
         % Make a little more room for the legends
         allylims(2)=diff(allylims)*1.5+allylims(1);
         for knd=1:6
           axes(ah(np+knd))
           % Plot reference diagonals as light gray lines on bottom
           rnd=0;
           refl=sqrt(10)./sqrt(2.^(2*(gv)));
           for mnd=allylims(1)-2:5:allylims(end)+45
             rnd=rnd+1;
             refd(rnd)=plot(gv,mnd+10*log10(refl),'-',...
                  'Color',[0.85 0.85 0.85],'HandleVisibility','off');
           end
           ah(np+knd).YLim=allylims;
         end
       end
       ti(ind)=title(sprintf('$%s$',labs{jnd}),'Interpreter','latex','Units','pixels');
       longticks
       xlim(minmax([glens glens2]).*[0.99 1.01])
       [axx(ind),xl(ind),yl(ind)]=...
         xtraxis(ah(ind),xxtks,xxtklbs,[]);
       axx(ind).XMinorTick='off';
       longticks
       ah(ind).XRuler.TickLabelGapOffset=-2;
       ah(ind).XTickLabelRotation=0;
       axx(ind).XRuler.TickLabelGapOffset=-3;
    else
       % Fourth row. Plot the rmse.
       rnd=ind-9;
       % The first dataset (100% observed grid)
       plot(glens,10*log10(abs(rmseemp(rnd,:))),'o',...
            'MarkerEdgeColor',clrs{1},'MarkerFaceColor','None',...
            'DisplayName',sprintf('100%% obs %0.1f',abs(mbfrmse(rnd,1))));
       plot(gv,predrmse(rnd,:),...
           'Color',clrs{1},'LineStyle','-','HandleVisibility','off');
       % The second dataset (66.7% observed grid)
       plot(glens2,10*log10(abs(rmseemp2(rnd,:))),'o',...
            'MarkerEdgeColor',clrs2{3},'MarkerFaceColor',clrs2{3},...
            'DisplayName',sprintf('66.7%% obs %0.1f',abs(mbfrmse2(rnd,1))));
       plot(gv2,predrmse2(rnd,:),...
            'Color',clrs2{3},'LineStyle','-','HandleVisibility','off');
       if rnd==np
         % Set consistent axes
         allylims=[ylim(ah([3*np+1:3*np+3]))];
         allylims=minmax([allylims{:}]);
         % Make a little more room for the legends
         allylims(2)=diff(allylims)*1.5+allylims(1);
         for knd=1:np
           axes(ah(3*np+knd))
           % Plot reference diagonals as light gray lines on bottom
           dnd=0;
           refl=sqrt(10)./sqrt(2.^(2*(gv)));
           for mnd=allylims(1)-2:5:allylims(end)+45
             dnd=dnd+1;
             refd(dnd)=plot(gv,mnd+10*log10(refl),'-',...
                  'Color',[0.85 0.85 0.85],'HandleVisibility','off');
           end
           ah(3*np+knd).YLim=allylims;
           if knd>1
             set(gca,'YTickLabel','')
           else
              ylabel(sprintf('$%s%s$',...
                     '10\times \mathrm{log}_{10}',...
                     '|\mathrm{rmse}\{\theta_i\}|'),...
                     'Interpreter','latex')
           end
           [leg1(3*np+knd),legic]=legend('interpreter','tex','EdgeColor','none',...
                                'NumColumns',2,'FontSize',7.5,...
                                'BackgroundAlpha',0,'box','off');
           % Have the legend markers take up less space
           movelcon(legic)
         end
       end
       axes(ah(ind))
       xlabel(xla,'Interpreter','latex','FontSize',fs)
       xlim(minmax([glens glens2]).*[0.99 1.01])
       ti(ind)=title(sprintf('$%s$',labs{ind}),'Interpreter','latex','Units','pixels');
       longticks
       xlim(minmax([glens glens2]).*[0.99 1.01])
       [axx(ind),xl(ind),yl(ind)]=...
         xtraxis(ah(ind),xxtks,xxtklbs,[]);
       axx(ind).XMinorTick='off';
       axx(ind).XRuler.TickLabelGapOffset=-3;
       longticks
     end
     axes(ah(ind))
     set(ah(ind),...
         'XTick',xtks,'XTickLabel',xtklbs,'XMinorTick','off',...
         'YMinorTick','off')
   end      
    
   % Super title to display the parameters
   % Include the geographic lengths of the grid in units of pi rho
   sti=sgtitle(...
              sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s$%0.1f%s%0.1f$','\sigma^2=',th1(1),...
              ' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),'m  $|$  dy=',...
              params.dydx(1),', dx=',params.dydx(2)),...
              'Interpreter','latex');

   % Shift the figure objects to make space for the title
   movev([ah(1:3) axx(1:3)],-0.03)
   movev([ah(4:12) axx(4:12)],-0.03)
   movev([ah([1:12]) axx([1:12])],-0.01)
   shrink([ah([1:12]) axx([1:12])],0.85,1)
   movev([leg1([1:12])],0.01)
   moveh(leg1([1:12]),0.01)
   movev(ti,12)
   movev(xl,0.1)
   
   % Final tinkering
   keyboard
   for ind=1:12
     if ind<10
       ah(ind).XTickLabel='';
     end
     axes(axx(ind))
     if ind>3
       axx(ind).XTickLabel='';
     end
   end
   movev([ah([4:12]) axx([4:12]) leg1([4:12])],0.015)
   movev([ah([7:12]) axx([7:12]) leg1([7:12])],0.015)
   movev([ah([10:12]) axx([10:12]) leg1([10:12])],0.015)
   moveh([ha([5:12]) axx([2:3,5:6,8:9,11:12]) leg1([2:3,5:6,8:9,11:12])],0.015)
   moveh([ha([9:12]) axx([3,6,9,12]) leg1([3,6,9,12])],0.015)
   for ind=1:12;ti(ind).FontSize=11;end

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
   np=3;

   % Labels
   labs1={'$\sigma^2$','$\nu$','$\rho$'};
   labs={labs1{:},'\sigma^2,\nu','\sigma^2,\rho', '\nu,\rho',labs1{:},labs1{:}};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];
   % Font size for latex interpreted text
   fs=11;

   % Set the Matern parameters 
   rh=16;
   th1=[3 0.5 rh];
   % Make a square grid with fixed grid length, gy
   gy=256;
   gygx=[gy gy];
   % Calculate factors of gy to choose multiple instances of NyNx and dydx with
   % the same gygx
   facs=1:ceil(sqrt(gy));
   dy=facs(rem(gy,facs)==0);
   ny=sort(gy./dy,'descend'); 
   % We should only take stock from experiments that have a sufficiently large
   % number of samples to begin with, which we will impose as 3 times pi*rho
   % (still pretty small)
   dy=[dy(ny>=pi*rh)];
   %dy=[0.5 dy(ny>=pi*rh*3)];
   ny=sort(gy./dy,'descend'); 
   ny=fliplr(ny);
   dy=fliplr(dy);
   Ns=ny.^2;
   % Confirm that each set of ny and dy results in the fixed grid size
   diferm(ny.*dy-gy)
   % Create these vectors for naming consistency with later plotting routines
   glens=log2(sqrt(Ns));
   log2glens=sqrt(Ns);
   % Second set of vectors for the second dataset
   glens2=log2(sqrt(Ns*2/3));
   log2glens2=sqrt(Ns*2/3);

   % We need a suite of estimates; set msk to full for a  unit taper, or random 
   % for a 66.7% observed field with the specified mask
   % In case we already calculated the estimates, we will try loading them

   % for random case, inc evenly in x and y: '11-May-2025';
   %datum='11-May-2025'; msk='random'; 
   datum='24-Aug-2025'; msk='random';
   datum='21-Jan-2026';msk='full';
   %'full';
   % fullx case, inc only in x: 
   %datum='14-May-2025'; msk='fullx';
   % full case, inc evenly in x and y: 
   %datum='01-May-2025'; msk='full';%'30-Apr-2025';
  % datum='24-Aug-2025';msk='full';

   % Flag for paper figure that will show the field fully observed and with a
   % constant ratio of random deletions on the same plot
   flg=1;
   incdft=0;
   if flg==1
     % We will let the fully observed field be loaded as usual
     datum='01-May-2025'; msk='full';
     datum='24-Aug-2025'; msk='full';
     datum='21-Jan-2026';msk='full';
     % We will load the field with random deletions here and denote its
     % variables with a '2'
     % 11-May-2025
     datum2='24-Aug-2025'; msk2='random';
     datum2='20-Jan-2026';msk2='random';
     demo7vars2=load(sprintf('covgammiosl_demo7_calcs_%s_%s.mat',datum2,msk2));
     thhats2=demo7vars2.thhats;
   %  cov1t2=demo7vars2.cov1t;
     cov3t2=demo7vars2.cov3t;
     covempt2=demo7vars2.covempt;
     ny2=demo7vars2.ny;
     dy2=demo7vars2.dy;
     params=demo7vars2.params;
     th1=demo7vars2.th1;
    % Ns2=sort((gy./dy2).^2,'descend');
     Ns2=Ns;
     % Distinguish each dataset by color
     clrs=["k","k","k"];
     clrs2=["#D95F02","#D95F02","#D95F02"];
   end
   try
   keyboard
     demo7calcs=load(sprintf('covgammiosl_demo7_calcs_%s_%s.mat',datum,msk));
     thhats=demo7calcs.thhats;
     covempt=demo7calcs.covempt;
     % cov1t=demo7calcs.cov1t;
     % try
     %   cov2t=demo7calcs.cov2t;
     %   if isnan(cov2t); incdft=0; else; incdft=1; end
     % catch
     %   incdft=0;
     % end
     cov3t=demo7calcs.cov3t;
     ny=demo7calcs.ny;
     dy=demo7calcs.dy;
     params=demo7calcs.params;
     th1=demo7calcs.th1;
   %  Ns=sort((gy./dy).^2,'descend');
   catch
     for ind=1:size(dy,2)
       % Set spacing and number of samples
       params.NyNx=[ny(ind) ny(ind)];
       params.dydx=[dy(ind) dy(ind)];
       if strcmp(msk,'random')
         % recalculate the mask taper based on grid size; hold constant for each
         % realization at this size so that the analytic covariance calculation
         % is specific to the same mask
         mcl=0.667;
         params.mask=msk;
         [~,I]=muckit(randn(params.NyNx),params,mcl);
         tapersv{ind}=I;
         params.taper=tapersv{ind};
       end
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
       mdout=0;
       if mdout
         % Apply MADOUT to robustly remove outliers; if outliers are removed,
         % continue simulating until acquiring NUMREALS estimates
         thhatst{ind}=madout(thhats(:,:,ind),1);
         while size(thhatst{ind},1)<numreals
           for jnd=1:(numreals-size(thhats,1))
             Hx=simulosl(th1,simparams,0);
             thhat=NaN;
             while isnan(thhat)
               [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
             end
             thhatspls(jnd,:)=thhat.*scl;
           end
           % Recreate thhats(:,:,ind) with the recalculated estimates
           thhats(:,:,ind)=[thhatst{ind};thhatspls];
           thhatspls=[];
           thhatst{ind}=madout(thhats(:,:,ind),1);
         end
         thhats(:,:,ind)=thhatst{ind};
       end
     end
     % Calculate parameter covariance up to 4 different ways. We may choose not to
     % include the dftmtx method in favor of looking at larger grids (the
     % diagonals method will provide identical results) 
     for ind=1:size(dy,2)
       % Set spacing and number of samples
       params.NyNx=[ny(ind) ny(ind)];
       params.dydx=[dy(ind) dy(ind)];
       if strcmp(msk,'random')
         params.taper=tapersv{:,ind};
       end
       % Calculate the empirical covariance from the ensemble of estimates
       covemp=nancov(thhats(:,:,ind));
       % Calculate parameter covariance according to the 'sample' method
       % cov1c=covthosl(th1,params,1,ifinv);
       % if incdft
       %   % Calculate parameter covariance according to the 'dftmtx' method
       %   cov2c=covthosl(th1,params,2,ifinv);
       % end
       % Calculate the exact covariance with the diagonals method
       meth=3;
       disp(sprintf('Calculating COVTHOSL for meth %i and NyNx %i %i',...
           meth,params.NyNx))
       cov3c=covthosl(th1,params,meth,ifinv);
       % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
       covempt(ind,:)=trilos(covemp);
       % if flg~=1; cov1t(ind,:)  =trilos(cov1c); end
       % if incdft; cov2t(ind,:)  =trilos(cov2c); end
       cov3t(ind,:)  =trilos(cov3c);
     end
     % Consider saving your calculations by setting 'svcalcs'
     svcalcs=1;
     if svcalcs
       if isfield(params,'mask')
         save(sprintf('covgammiosl_demo7_calcs_%s_%s.mat',datum,msk),...
              'thhats','cov3t','covempt','ny','dy','params','th1','tapersv')
       else
         save(sprintf('covgammiosl_demo7_calcs_%s_%s.mat',datum,msk),...
              'thhats','cov3t','covempt','ny','dy','params','th1')
       end
     end
   end

   clf
   figip=get(gcf,'InnerPosition');
   set(gcf,'InnerPosition',[figip(1:2) 950 950])
   [ah,ha,H]=krijetem(subnum(4,3));
   fs=11;
   xla='$\sqrt{K}$';
   gv=0.8*glens(1):0.1:1.2*glens(end);
   gv2=0.8*glens2(1):0.1:1.2*glens2(end);
   % I used to be interested in showing that each of the covthosl calculations
   % are essentially identical, but since this is established, we won't bother
   % to muddy the plot:
   covana=cov3t;
   covana2=cov3t2;
   % Complete the calculations. Fit a slope to the covariance calculations as a 
   % function of grid length 
   for jnd=1:6
     mbfemp(jnd,:)     =polyfit(glens, 10*log10(abs(covempt(:,jnd))),1);
     mbfcovana(jnd,:)  =polyfit(glens, 10*log10(abs(covana(:,jnd))),1);
     predemp(jnd,:)    =polyval(mbfemp(jnd,:),   gv);
     predcovana(jnd,:) =polyval(mbfcovana(jnd,:),  gv);
     mbfemp2(jnd,:)    =polyfit(glens2,10*log10(abs(covempt2(:,jnd))),1);
     mbfcovana2(jnd,:) =polyfit(glens2,10*log10(abs(covana2(:,jnd))),1);
     predemp2(jnd,:)   =polyval(mbfemp2(jnd,:),   gv2);
     predcovana2(jnd,:)=polyval(mbfcovana2(jnd,:),gv2);
     if jnd<4
       % rmse 
       rmseemp(jnd,:)   =sqrt(squeeze(var(thhats(:,jnd,:)))+...
                             squeeze(mean(thhats(:,jnd,:))-th1(jnd)).^2);
       mbfrmse(jnd,:)   =polyfit(glens,10*log10(abs(rmseemp(jnd,:))),1);
       predrmse(jnd,:)  =polyval(mbfrmse(jnd,:),   gv);
       rmseemp2(jnd,:)  =sqrt(squeeze(var(thhats2(:,jnd,:)))+...
                             squeeze(mean(thhats2(:,jnd,:))-th1(jnd)).^2);
       mbfrmse2(jnd,:)  =polyfit(glens2,10*log10(abs(rmseemp2(jnd,:))),1);
       predrmse2(jnd,:) =polyval(mbfrmse2(jnd,:),   gv2);
     end
   end

   for ind=1:12
     % First row. Plot the relative bias of the estimates from the empirical study
     axes(ah(ind));
     hold on
     if ind<4
       relbsme(:,ind)=squeeze(mean(thhats(:,ind,:)-th1(ind)))/th1(ind);
       relbssd(:,ind)=squeeze(sqrt(var(thhats(:,ind,:)-th1(ind),"omitnan"))/...
                      th1(ind));
       relbsme2(:,ind)=squeeze(mean(thhats2(:,ind,:)-th1(ind)))/th1(ind);
       relbssd2(:,ind)=squeeze(sqrt(var(thhats2(:,ind,:)-th1(ind),"omitnan"))/...
                      th1(ind));
       % First dataset (100% observed grid)
       % Errorbar in terms of the standard deviation 
       pl1(ind)=errorbar(glens,relbsme(:,ind),relbssd(:,ind),'ko','CapSize',5,...
                   'DisplayName','100% obs');
       % Add the bias markers for the second dataset
       pl12(ind)=errorbar(glens2,relbsme2(:,ind),relbssd2(:,ind),'o',...
                   'MarkerEdgeColor',clrs2{1},'MarkerFaceColor',clrs2{1},...
                   'LineWidth',1,'CapSize',5,'Color',clrs2{1},...
                   'DisplayName','66.7% obs');
       xlim(minmax([glens glens2]).*[0.99 1.01])
      % We don't want this anymore
      % if ind==2
      %   % the bias in nu for the full grid is much smaller, so lets pull it on
      %   % top
      %   uistack(pl1(2),'up',1)
      % end
       % Add the 0-line where the estimate equals the truth
       yline(0,'HandleVisibility','off')
       if ind==1
         ylabel('$(\langle\hat{\theta}\rangle-\theta_0)/\theta_0$','interpreter','latex');
       elseif ind==2
         ylim([-0.1,max(nanmean(thhats(:,ind,:))-th1(ind))+0.1])
       end 
       ti(ind)=title(labs{ind},'Interpreter','latex','Units','pixels');
       [leg1(ind),legic]=legend('interpreter','tex','EdgeColor','none','box','on',...
                            'NumColumns',2,'FontSize',7.5,...
                            'BackgroundAlpha',0,'box','off');
       % Have the legend markers take up less space
       movelcon(legic)
       
       if ind==3
         % Cosmetics for the first row, including second xaxis
         byl=max(max(abs(cell2mat(ylim(ah([1:3]))))));
         byl=[-byl byl];
         byl=[-0.8 0.8];
         for knd=1:ind
           axes(ah(knd))
           ylim([byl])
           box on
           longticks
           if knd>1
             set(gca,'YTickLabel','')
           end
           % Set the top x-axis
           xlb='$\Delta y = \Delta x$ Spacing';
           [axx(knd),xl(knd),yl(knd)]=...
             xtraxis(ah(knd),xtks,dy,xlb);
           set(xl(knd),'Interpreter','latex')
           longticks
           ah(knd).XRuler.TickLabelGapOffset=-2;
           ah(knd).XTickLabelRotation=0;
           axx(knd).XRuler.TickLabelGapOffset=-3;
         end
       end
     elseif ind<10
       % The first and second rows. Plot the absolute (co)variance.
       % The first dataset (100% observed grid)
       jnd=ind-np;
       pe{1,jnd}=plot(glens,10*log10(abs(covempt(:,jnd))),'ko',...
                       'DisplayName',sprintf('100%% obs %0.1f',abs(mbfemp(jnd,1))));
       pa{1,jnd}=plot(glens,10*log10(abs(covana(:,jnd))) ,'^','MarkerEdgeColor',clrs{1},...
              'MarkerFaceColor','None',...
              'DisplayName',sprintf('100%% pred %0.1f',abs(mbfcovana(jnd,1))));
       % Plot the absolute covariance calculated from the ensemble and covthosl
       ppe{1,jnd}=plot(gv,predemp(jnd,:),'Color','k','HandleVisibility','off');
       ppa{2,jnd}=plot(gv,predcovana(jnd,:),'Color',clrs{1},'HandleVisibility','off');

       % The second dataset (66.7% observed grid, random deletions)
       pe{2,jnd}=plot(glens2,10*log10(abs(covempt2(:,jnd))),'o','MarkerEdgeColor',clrs2{1},...
                      'MarkerFaceColor',clrs2{1},...
                      'DisplayName',sprintf('66.7%% obs %0.1f',abs(mbfemp2(jnd,1))));
       pa{2,jnd}=plot(glens2,10*log10(abs(covana2(:,jnd))),'^','MarkerEdgeColor',clrs2{1},...
                 'MarkerFaceColor',clrs2{1},...
                 'DisplayName',sprintf('66.7%% pred %0.1f',abs(mbfcovana2(jnd,1))));
       % Plot the absolute covariance calculated from the ensemble and covthosl
       ppe{1,jnd}=plot(gv2,predemp2(jnd,:),'Color',clrs2{3},'LineStyle','-',...
                       'HandleVisibility','off');
       ppa{2,jnd}=plot(gv2,predcovana2(jnd,:),'Color',clrs2{3},'LineStyle','-',...
                       'HandleVisibility','off');
       %set(gca,'YLim',[1e-5 1e7],'YTick',[1e-3 1e0 1e3 1e6])
       if ~ismember(jnd,[1 4])
         set(gca,'YTickLabel','')
       end
       [leg1(ind),legic]=legend('interpreter','tex','EdgeColor','none','box','on',...
                                'NumColumns',2,'FontSize',7.5,...
                                'BackgroundAlpha',0,'box','off');
       leg1(ind).AutoUpdate='off';
       % Have the legend markers take up less space
       movelcon(legic)
       
       % Labels, titles, ticks, common limits
       if jnd==1
          ylabel('$10\times \mathrm{log}_{10}(|\mathrm{var}\{\theta_i\}|)$',...
               'Interpreter','latex','FontSize',fs)
       elseif jnd==4
          ylabel('$10\times \mathrm{log}_{10}(|\mathrm{cov}\{\theta_i,\theta_j\}|$)',...
                 'Interpreter','latex','FontSize',fs)
       end
       if jnd==6
         % Set consistent axes
         allylims=[ylim(ah([np+1:np+6]))];
         allylims=minmax([allylims{:}]);
         % Make a little more room for the legends
         allylims(2)=diff(allylims)*1.5+allylims(1);
         for knd=1:6
           axes(ah(np+knd))
           % Plot reference diagonals as light gray lines on bottom
           rnd=0;
           refl=sqrt(10)./sqrt(2.^(2*(gv)));
           for mnd=allylims(1)-2:5:allylims(end)+45
             rnd=rnd+1;
             refd(rnd)=plot(gv,mnd+10*log10(refl),'-',...
                  'Color',[0.85 0.85 0.85],'HandleVisibility','off');
           end
           ah(np+knd).YLim=allylims;
         end
       end
       ti(ind)=title(sprintf('$%s$',labs{jnd}),'Interpreter','latex','Units','pixels');
       longticks
       box on
       xlim(minmax([glens glens2]).*[0.99 1.01])
       [axx(ind),xl(ind),yl(ind)]=...
         xtraxis(ah(ind),xtks,dy,[]);
       longticks
       box on
       ah(ind).XRuler.TickLabelGapOffset=-2;
       ah(ind).XTickLabelRotation=0;
       if gy==120
         ah(ind).XTickLabel(3,:)=' ';
       end
       axx(ind).XRuler.TickLabelGapOffset=-3;
     else
       % Fourth row. Plot the rmse.
       rnd=ind-9;
       % The first dataset (100% observed grid)
       plot(glens,10*log10(abs(rmseemp(rnd,:))),'o',...
            'MarkerEdgeColor',clrs{1},'MarkerFaceColor','None',...
            'DisplayName',sprintf('100%% obs %0.1f',abs(mbfrmse(rnd,1))));
       plot(gv,predrmse(rnd,:),...
           'Color',clrs{1},'LineStyle','-','HandleVisibility','off');
       % The second dataset (66.7% observed grid)
       plot(glens2,10*log10(abs(rmseemp2(rnd,:))),'o',...
            'MarkerEdgeColor',clrs2{3},'MarkerFaceColor',clrs2{3},...
            'DisplayName',sprintf('66.7%% obs %0.1f',abs(mbfrmse2(rnd,1))));
       plot(gv2,predrmse2(rnd,:),...
            'Color',clrs2{3},'LineStyle','-','HandleVisibility','off');
       if rnd==np
         % Set consistent axes
         allylims=[ylim(ah([3*np+1:3*np+3]))];
         allylims=minmax([allylims{:}]);
         % Make a little more room for the legends
         allylims(2)=diff(allylims)*1.5+allylims(1);
         for knd=1:np
           axes(ah(3*np+knd))
           % Plot reference diagonals as light gray lines on bottom
           dnd=0;
           refl=sqrt(10)./sqrt(2.^(2*(gv)));
           for mnd=allylims(1)-2:5:allylims(end)+45
             dnd=dnd+1;
             refd(dnd)=plot(gv,mnd+10*log10(refl),'-',...
                  'Color',[0.85 0.85 0.85],'HandleVisibility','off');
           end
           ah(3*np+knd).YLim=allylims;
           if knd>1
             set(gca,'YTickLabel','')
           else
              ylabel(sprintf('$%s%s$',...
                     '10\times \mathrm{log}_{10}',...
                     '|\mathrm{rmse}\{\theta_i\}|'),...
                     'Interpreter','latex')
           end
           [leg1(3*np+knd),legic]=legend('interpreter','tex','EdgeColor','none','box','on',...
                                'NumColumns',2,'FontSize',7.5,...
                                'BackgroundAlpha',0,'box','off');
           % Have the legend markers take up less space
           movelcon(legic)
         end
       end
       axes(ah(ind))
       xlabel(xla,'Interpreter','latex','FontSize',fs)
       xlim(minmax([glens glens2]).*[0.99 1.01])
       ti(ind)=title(sprintf('$%s$',labs{ind}),'Interpreter','latex','Units','pixels');
       longticks
       box on
       xlim(minmax([glens glens2]).*[0.99 1.01])
       xtks=glens;
       xtklbs=string(log2glens);
       set(ah(ind),...
           'XTick',xtks,'XTickLabel',xtklbs,'XMinorTick','off',...
           'YMinorTick','off')
       [axx(ind),xl(ind),yl(ind)]=...
         xtraxis(ah(ind),xtks,dy,[]);
       axx(ind).XRuler.TickLabelGapOffset=-3;
       axx(ind).XRuler.TickLabelGapOffset=-3;
       longticks
     end
     xtks=glens;
     xtklbs=string(log2glens);
     set(ah(ind),...
         'XTick',xtks,'XTickLabel',xtklbs,'XMinorTick','off',...
         'YMinorTick','off','Box','on')
   end      
    
   % Super title to display the parameters
   % Include the geographic lengths of the grid in units of pi rho
   gypr=round(gy/pi/th1(3));
   sti=sgtitle(...
              sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s$%i%s%i%s$','\sigma^2=',th1(1),...
              ' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),'m  $|$  ',...
              gypr,'[\pi \rho] \times',gypr,'[\pi \rho]'),...
              'Interpreter','latex');
   % Shift the figure objects to make space for the title
   movev([ah(1:3) axx(1:3)],-0.03)
   movev([ah(4:12) axx(4:12)],-0.03)
   movev([ah([1:12]) axx([1:12])],-0.01)
   shrink([ah([1:12]) axx([1:12])],0.85,1)
   movev([leg1([1:12])],0.01)
   moveh(leg1([1:12]),0.01)
   movev(ti,12)
   movev(xl,0.1)

   % Final tinkering
   for ind=1:12
     if ind<10
       ah(ind).XTickLabel='';
     end
     axes(axx(ind))
     if ind>3
       axx(ind).XTickLabel='';
     end
   end
   movev([ah([4:12]) axx([4:12]) leg1([4:12])],0.015)
   movev([ah([7:12]) axx([7:12]) leg1([7:12])],0.015)
   movev([ah([10:12]) axx([10:12]) leg1([10:12])],0.015)
   moveh([ha([5:12]) axx([2:3,5:6,8:9,11:12]) leg1([2:3,5:6,8:9,11:12])],0.015)
   moveh([ha([9:12]) axx([3,6,9,12]) leg1([3,6,9,12])],0.015)
   for ind=1:12;ti(ind).FontSize=11;end
   for ind=4:9;ah(ind).YLim=[-40 60];end

   keyboard
   % Save the figure
   saveas(gcf,...
     sprintf('covgammiosl_demo7_infill_%inu%irh%i_%s.eps',...
       th1.*[1 10 1],date),'epsc')

 %  % Make a second figure to show w(x) and fft2(|w(x)|^2) for changes in 
 %  % spacings for the random deletion case
 %  figure();
 %  clf
 %  set(gcf,'InnerPosition',[1 1 950 200])
 %  [ah,ha]=krijetem(subnum(1,6));
 % % % Plot an example of the 66.7% taper for the first grid
 % % axes(ah(1)); 
 % % imagesc(demo7vars2.tapersv{1})
 % % axis image; longticks

 % % % Plot the kernel of the 66.7% taper for the first grid
 % % axes(ah(2));
 % % imagesc(log10(abs(fft2(demo7vars2.tapersv{1})).^2))
 % % axis image; longticks

 % % % Plot an example of the 66.7% taper for the first grid
 % % axes(ah(3));
 % % imagesc(demo7vars2.tapersv{2})
 % % axis image; longticks

 % % % Plot the kernel of the 66.7% taper for the first grid
 % % axes(ah(4));
 % % imagesc(log10(abs(fft2(demo7vars2.tapersv{2})).^2))
 % % axis image; longticks

 % % % Plot an example of the 66.7% taper for the first grid
 % % axes(ah(5));
 % % imagesc(demo7vars2.tapersv{3})
 % % axis image; longticks

 % % % Plot the kernel of the 66.7% taper for the first grid
 % % axes(ah(6));
 % % imagesc(log10(abs(fft2(demo7vars2.tapersv{3})).^2))
 % % axis image; longticks

 %  %ftp=abs(fft2(demo7vars2.tapersv)).^2;
 %  % Let's take a quick look at each taper 
 %  for ind=1:numel(dy)
 %    ftx0=fftshift(fft2(demo7vars2.tapersv{ind}));
 %    ftph0=imag(ftx0);
 %    ftph{ind}=ftph0;
 %    ftx0=decibel(abs(ftx0).^2);
 %    % Decibel will produce -Inf if an element is 0; replace these with the
 %    % next smallest value
 %    ftx0(isinf(ftx0))=NaN;
 %    ftx0(isnan(ftx0))=min(ftx0,[],'all');
 %    ftx{ind}=ftx0;
 %    % And again, but this time with padding
 %    ftx0p=fftshift(fft2(demo7vars2.tapersv{ind},2^nextpow2(ny(ind)),...
 %                                                2^nextpow2(ny(ind))));
 %    ftph0p=imag(ftx0p);
 %    ftphp{ind}=ftph0p;
 %    ftx0p=decibel(abs(ftx0p).^2);
 %    % Decibel will produce -Inf if an element is 0; replace these with the
 %    % next smallest value
 %    ftx0p(isinf(ftx0p))=NaN;
 %    ftx0p(isnan(ftx0p))=min(ftx0p,[],'all');
 %    ftxp{ind}=ftx0p;
 %  end
 %  % Plot with padding?
 %  upad=0;
 %  if upad; ftp=ftxp; ftphs=ftphp; else; ftp=ftx; ftphs=ftph; end

 %  axs=1:size(ha,2);
 %  xnds=axs(~~mod(axs,2));
 %  knds=axs(~mod(axs,2));
 %  for ind=1:size(ha,2)/2
 %    % Display the taper in grid space
 %    xnd=xnds(ind);
 %    knd=knds(ind);
 %    axes(ha(xnd))
 %    imagesc(demo7vars2.tapersv{ind})
 %    colormap(gca,gray(2))
 %    caxis([0 1])
 %    hold on
 %    % Show the decorrelation length scale (the correlation is 1/3 of the
 %    % variance at a distance of pi*rho, unitless)
 %    dc(ind,1:2)=circ(pi*th1(3)/dy(ind),[],[pi*th1(3)/dy(ind) ny(ind)-pi*th1(3)/dy(ind)]); 
 %    set(dc(ind,1),'color','b'); delete(dc(ind,2))
 %    ind=1
 %    tti(ind)=title(sprintf('$%sx=%i, Nx=%i:%s$','$\Delta$',...
 %                           dy(ind),ny(ind)),'$w(\mathbf{x})$',...
 %                   'Interpreter','latex');
 %    longticks
 %    axis image
 %    shrink(ha(xnd),1.05,1.05)
 %    set(ha(xnd),'XTick',[1 ny(ind)],'YTick',[1 ny(ind)])
 %    xlabel('$x_1$','Interpreter','latex')
 %    ylabel('$x_2$','Interpreter','latex')

 %    % Display the taper in wavenumber space
 %    axes(ha(knd))
 %    colormap(gca,parula)
 %    pnd=1;
 %    if pnd==1
 %      % square modulus pf the FT of the taper
 %      imagesc(ftp{ind})
 %      caxis(minmax([ftp{1}(:); ftp{2}(:); ftp{3}(:)]))
 %      tti2(ind)=title('$|\mathrm{fft}(w(\mathbf{x})|^2$','Interpreter','latex');
 %    else
 %      % Phase of the FT of the taper
 %      imagesc(ftphs{ind})
 %      %caxis(minmax(ftphs(:,:,2:end)))
 %    end
 %    longticks
 %    axis image
 %    shrink(ha(knd),1.05,1.05)
 %    if upad==1
 %      set(ha(knd),'XTick',[1 ny(ind)],'YTick',[1 ny(ind)])
 %    else
 %      set(ha(knd),'XTick',[1 size(ftxp{ind},1)],'YTick',[1 size(ftxp{ind},2)])
 %    end
 %    xlabel('$k_1$','Interpreter','latex')
 %    ylabel('$k_2$','Interpreter','latex')
 %    set(ha(knd),'YAxisLocation','right')

 %    moveh([ha(xnd) ha(knd)],-0.01*(3-4*(ind-1)))
 %    moveh(ha(knd),-0.02)
 %  end
elseif strcmp(th,'demo8')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the PERCENT OF MISSINGNESS, either due to
   % random deletions or a growing mask

   % Option to display multiple experiments at once that use the
   % same mask, specifically for increasing values of rho. Color code connecting
   % line segment with value of rho; use same color in decorrelation radius
   % circle!

   % The Matern and regular grid parameters
   th1=[3 0.5 10];
   params.NyNx=[88 88];params.dydx=[1 1];
   %params.NyNx=[128 128];params.dydx=[1 1];
   params.blurs=-1;params.taper=1;                                            
   ifinv=[1 1 1];                                                              
   oneprcnt=ceil(0.01.*prod(params.NyNx));
   % If you specify a mask, it needs to be usable by MASKIT, either as an index
   % matrix or a region name that we have worked with previously
   %params.mask='puss';
   %params.mask='random';datum='30-Apr-2025';% for random
   %params.mask='random';datum='30-Nov-2025';% for random
   % CHANGED ON 12/23 
   %params.mask='random';datum='03-Dec-2025';% for random
   params.mask='random';datum='23-Dec-2025';% for random
   datum='20-Jan-2026';
   % remove
   %params.mask='france';datum='30-Apr-2025'; 
   % france for rh=10: '15-May-2025';
   % france for rh=8: '14-May-2025';
   % france for rh=6: '13-May-2025';
   % france for rh=4: '09-May-2025';%'30-Apr-2025'
   %params.mask='Italy';datum='08-May-2025';
   %params.mask='doves';datum='11-May-2025';
   %params.mask='newjersey';datum='11-May-2025';
   if isfield(params,'mask')
     % If you provided a mask, randomly select the grid index location that the
     % mask will be centered on (p.mymx) to be passed onto MASKIT
     vrshft=randn*0.2*params.NyNx(1)/2;
     hzshft=randn*0.2*params.NyNx(2)/2;
     vrshft=0;
     hzshft=0;
     params.mymx=[vrshft hzshft];
   end
   % Simulate via circulant embedding
   simparams=params;simparams.blurs=Inf;
   % If you did not specify a mask, we will proceed with random deletions
   % The number of parameters
   np=3;

   glens=3.9:0.5:6.4;
   log2glens=2.^glens;
   Ks=ceil(log2glens);
   prcobs=floor(Ks./sqrt(prod(params.NyNx)).*100);
   missarr=(100-prcobs)./100;
   % Make two different plot styles; if flg==1, publication style figure
   flg=1;
   if flg==1
     th1=[3 0.5 2];
     ds1.th1=th1;ds1.th1(3)=2;
     ds2.th1=th1;ds2.th1(3)=4;
     ds3.th1=th1;ds3.th1(3)=8;
     try
       ds1=load(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',ds1.th1(3),datum));
       ds2=load(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',ds2.th1(3),datum));
       ds3=load(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',ds3.th1(3),datum));
       numreals=size(ds1.thhats,1);
     catch
       numreals=100;
       % random missingness
       thhats1=zeros(numreals,np,numel(missarr));
       thhats2=zeros(numreals,np,numel(missarr));
       thhats3=zeros(numreals,np,numel(missarr));
       taper=ones(params.NyNx);
       lstmissarr=0;
       mktprmarr=fliplr(missarr);
       for ind=1:numel(mktprmarr)
         while sum(taper(:))./prod(params.NyNx)>1-mktprmarr(ind)
           taper(randi(prod(params.NyNx),1,...                                
                       ceil(oneprcnt.*(abs(1-mktprmarr(ind)-lstmissarr)))))=deal(0);
                       %100.*ceil(oneprcnt.*(abs(missarr(ind)-lstmissarr)))))=deal(0);
           lstmissarr=sum(taper(:))./prod(params.NyNx);
         end
         tapersv(:,:,numel(mktprmarr)-ind+1)=taper;
       end
       for ind=1:numel(missarr)
         params.taper=tapersv(:,:,ind);simparams.taper=tapersv(:,:,ind);
         for mnd=1:numreals
           thhat=NaN;
           while isnan(thhat)                                                   
             Hx=simulosl(ds1.th1,simparams,0);                                      
             [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],ds1.th1,ifinv,0);            
           end                                                                  
           thhats1(mnd,:,ind)=deal(thhat.*scl);                                  

           thhat=NaN;
           while isnan(thhat)                                                   
             Hx=simulosl(ds2.th1,simparams,0);                                      
             [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],ds2.th1,ifinv,0);            
           end                                                                  
           thhats2(mnd,:,ind)=deal(thhat.*scl);                                  

           thhat=NaN;
           while isnan(thhat)                                                   
             Hx=simulosl(ds3.th1,simparams,0);                                      
             [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],ds3.th1,ifinv,0);            
           end                                                                  
           thhats3(mnd,:,ind)=deal(thhat.*scl);                                  
         end                                                                    
       end
       keyboard
       for ind=1:numel(missarr)                                                 
         params.taper=tapersv(:,:,ind);                                         
         covemp1=nancov(thhats1(:,:,ind));                                        
         covemp2=nancov(thhats2(:,:,ind));                                        
         covemp3=nancov(thhats3(:,:,ind));                                        
         try
           meth=2;
           disp(sprintf('Calculating COVTHOSL for meth %i and NyNx %i %i',...
             meth,params.NyNx))
           covc1=covthosl(ds1.th1,params,meth,ifinv);
           covc2=covthosl(ds2.th1,params,meth,ifinv);
           covc3=covthosl(ds3.th1,params,meth,ifinv);
         catch
          meth=3;
          disp(sprintf('Calculating COVTHOSL for meth %i and NyNx %i %i',...
            meth,params.NyNx))
          covc1=covthosl(ds1.th1,params,meth,ifinv);
          covc2=covthosl(ds2.th1,params,meth,ifinv);
          covc3=covthosl(ds3.th1,params,meth,ifinv);
         end

         cov3t1(ind,:)=trilos(covc1);                                            
         covet1(ind,:)=trilos(covemp1);                                           

         cov3t2(ind,:)=trilos(covc2);                                            
         covet2(ind,:)=trilos(covemp2);                                           

         cov3t3(ind,:)=trilos(covc3);                                            
         covet3(ind,:)=trilos(covemp3);                                           
       end                                                                      
       keyboard
       cov3t=cov3t1;
       covet=covet1;
       thhats=thhats1;
       th1=ds1.th1;
       %% TODO: save the additional estimates to the existing file (12/1)!
       save(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',ds1.th1(3),datum),...   
            'thhats','tapersv','cov3t','covet','missarr','th1')

       cov3t=cov3t2;
       covet=covet2;
       thhats=thhats2;
       th1=ds2.th1;
       save(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',ds2.th1(3),datum),...   
            'thhats','tapersv','cov3t','covet','missarr','th1')

       cov3t=cov3t3;
       covet=covet3;
       thhats=thhats3;
       th1=ds3.th1;
       save(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',ds3.th1(3),datum),...   
            'thhats','tapersv','cov3t','covet','missarr','th1')
     end
   else
     try
       % Load experiments
       if isfield(params,'mask')
         keyboard
         load(sprintf('covgammiosl_demo8vars_%s_%s.mat',params.mask,datum),...
              'thhats','tapersv','cov1t','cov2t','cov3t','covet','missarr','prcmiss')
         if all(isnan(cov2t)); incdft=0; else; incdft=1; end
         numreals=size(thhats,1);
       else
         load(sprintf('covgammiosl_demo8vars_random_%s.mat',datum),...
              'thhats','tapersv','cov1t','cov2t','cov3t','covet','missarr','prcmiss')
         if all(isnan(cov2t)); incdft=0; else; incdft=1; end
         numreals=size(thhats,1);
       end
     catch
       % How many realizations for calculating the empirical covariance?
       numreals=200;
       % What is the maximum percent missingness we will study? If we are working
       % with masks, we need to remember that the scale we provide is not directly
       % proportional to the amount of missing observations
       if isfield(params,'mask')
         if strcmp(params.mask,'france')
           % max missing will be 56.5%
           missarr=logspace(0,2,9);
         elseif strcmp(params.mask,'Italy')
           % max missing will be 56.55%
           missarr=logspace(0,2.85,9);
         elseif strcmp(params.mask,'doves')
           missarr=logspace(0,2.196,9);
         elseif strcmp(params.mask,'newjersey')
           missarr=logspace(0,2.057,9);
         end
       else
         % random deletions
         missarr=logspace(0,2,9);
       end
       % Preallocate space for the estimates
       thhats=zeros(numreals,np,numel(missarr));
       % Calculate the estimates for increasing missingness for the mask type
       % selected 
       taper=ones(params.NyNx);
       for ind=1:numel(missarr)
         params.taper=taper;simparams.taper=taper;
         tapersv(:,:,ind)=taper;
         % Calculate the percentage of the grid that is missing
         prcmiss(ind)=1-sum(taper,"all")./prod(params.NyNx);
         prcobs(ind)=sum(taper,"all")./prod(params.NyNx);
         for mnd=1:numreals
           thhat=NaN;
           while isnan(thhat)
             Hx=simulosl(th1,simparams,0);
             [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
           end
           thhats(mnd,:,ind)=deal(thhat.*scl);
         end
         if ind<numel(missarr)
           if isfield(params,'mask')
             [~,I]=maskit(rand(params.NyNx),params,missarr(ind+1)/100);
             taper=~I;
             if strcmp(params.mask,'doves')
               taper=flipud(fliplr(taper));
             end
           else
             taper(randi(prod(params.NyNx),1,...
                         ceil(oneprcnt.*(missarr(ind+1)-missarr(ind)))))=deal(0);
                         %ceil((missarr(ind+1)-missarr(ind))/100)))=deal(0);
           end
         end
       end
       % Calculate parameter covariance up to 4 different ways and plot as a
       % function of the percentage of the regular grid that is missing. We may 
       % choose not to include the dftmtx method so that we can look at larger grids
       incdft=0;
       for ind=1:numel(missarr)
         params.taper=tapersv(:,:,ind);
         covemp=nancov(thhats(:,:,ind));
         covg1 =covgammiosl(th1,params,1,ifinv);
         if incdft; covg2 =covgammiosl(th1,params,2,ifinv); end
         covg3 =covgammiosl(th1,params,3,ifinv);
         cov1c=varianceofestimates(th1,params,covg1,ifinv);
         if incdft; cov2c=varianceofestimates(th1,params,covg2,ifinv); end
         cov3c=varianceofestimates(th1,params,covg3,ifinv);
         cov1t(ind,:)=trilos(cov1c);
         if incdft; cov2t(ind,:)=trilos(cov2c); else; cov2t(ind,:)=NaN; end
         cov3t(ind,:)=trilos(cov3c);
         covet(ind,:)=trilos(covemp);
       end
       save(sprintf('covgammiosl_demo8vars_%s_%s.mat',params.mask,date),...
            'thhats','tapersv','cov1t','cov2t','cov3t','covet','missarr','prcmiss','th1')
     end
   end

   mdrmv=0;
   if mdrmv
     % Apply MADOUT to robustly remove outliers; if outliers are removed,
     % continue simulating until acquiring NUMREALS estimates
     for ind=1:numel(missarr)
       thhatst{ind}=madout(thhats(:,:,ind),1);
       params.taper=tapersv(:,:,ind);
       while size(thhatst{ind},1)<numreals
         for jnd=1:(numreals-size(thhatst{ind},1))
           Hx=simulosl(th1,simparams,0);
           thhat=NaN;
           while isnan(thhat)
             [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
           end
           thhatspls(jnd,:)=thhat.*scl;
         end
         % Recreate thhats(:,:,ind) with the recalculated estimates
         thhats(:,:,ind)=[thhatst{ind};thhatspls];
         thhatspls=[];
         thhatst{ind}=madout(thhats(:,:,ind),1);
       end
       thhats(:,:,ind)=thhatst{ind};
     end
     covemp=nancov(thhats(:,:,ind));
     covet(ind,:)=trilos(covemp);
   end

   % Begin the figure
   % Labels
   labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
         '\sigma^2,\rho', '\nu,\rho'};
   if ~flg
     % Color friendly palette
     clrs=["#1B9E77","#D95F02","#7570B3"];

     clf
     %[ah,ha,H]=krijetem(subnum(3,3));
     [ah,ha,H]=krijetem(subnum(4,3));

     % Plot the relative bias of the estimates from the empirical study
     labs2={'$\sigma^2$','$\nu$','$\rho$'};
     for ind=1:np
       rnd=ind;
       axes(ah(rnd));
       relbsme(:,ind)=squeeze(nanmean(thhats(:,ind,:)-th1(ind)))/th1(ind);
       relbssd(:,ind)=squeeze(sqrt(nanvar(thhats(:,ind,:)-th1(ind),"omitnan"))/...
                      th1(ind));
       % Errorbar in terms of the standard deviation 
       pl1(ind)=errorbar(prcmiss*100,relbsme(:,ind),relbssd(:,ind),'ko');
       pl1(ind).CapSize = 0;
       hold on
       % Add the 0-line where the estimate equals the truth
       yline(0)
       xlabel('Percent missing')
       if ind==1
         ylabel('$(\langle\hat{\theta}\rangle-\theta_0)/\theta_0$','interpreter','latex');
       end
       ti(rnd)=title(labs2{ind},'Interpreter','latex');
       longticks
      % xlim([min(prcmiss*100) 100]+[-0.5 0.5])
      xlim([0.9 100])
      set(gca,'XScale','log','XTick',[0 1 10 20 40 100])%,...
%              'XTickLabel',{'1';'10';'20';'40','100'})
     end

     %byl=max(abs([ylim(ah(7)) ylim(ah(8)) ylim(ah(9))]));
     byl=max(abs([ylim(ah(1)) ylim(ah(2)) ylim(ah(3))]));
     byl=[-byl byl];
     % for consistency between the france and random experiment
     byl=[-0.5 0.5];
     for ind=1:np
       %axes(ah(ind+6))
       axes(ah(ind))
       set(gca,'YLim',byl,'XMinorTick','off')          
     end

     if isfield(params,'mask')
       masknm=params.mask;
     else
       masknm='random';
     end
     for ind=1:numel(missarr)
       for jnd=1:6
         rnd=np+jnd;
         axes(ah(rnd))
         plot(prcmiss(ind)*100,abs(covet(ind,jnd)),'ko','DisplayName','empirical')
         hold on
         plot(prcmiss(ind)*100,abs(cov1t(ind,jnd)),'^','MarkerEdgeColor',clrs(1),...
              'MarkerFaceColor',clrs(1),'DisplayName','sample')
         if incdft; plot(prcmiss(ind)*100,abs(cov2t(ind,jnd)),'s','MarkerEdgeColor',clrs(2),...
              'MarkerFaceColor',clrs(2),'DisplayName','dftmtx'); end
         plot(prcmiss(ind)*100,abs(cov3t(ind,jnd)),'*','MarkerEdgeColor',clrs(3),...
              'MarkerFaceColor',clrs(3),'DisplayName','diagonals')
         if jnd==2 & ind==1
           % Make a legend the first time we touch the top center plot
           if incdft
             [leg1,legic]=legend('empirical','sample','dftmtx','diagonals',...
                                 'box','off');
           else
             [leg1,legic]=legend('empirical','sample','diagonals',...
                                 'box','off');
           end
           leg1.AutoUpdate='off';
           % Reduce the amount of space each legend marker uses 
           movelcon(legic)
           % Move the legend to the left edge of the axis
           leg1.Position(1)=0.39;
           leg1.Position(2)=0.83;
        end
      end
     end
     % for consistency between experiments that we are studying, set same limits
     ayl={[0.1 1.2],[0.000001 0.03],[0.4 2.5],[0.001 0.015],[0.1 2.1],[0.005 0.08]};
     for jnd=1:6
      axes(ah(jnd))
      if jnd==1|jnd==4
         ylabel('$|\mathrm{cov}\{\theta_i,\theta_j\}|$','Interpreter','latex')
      end
      title(sprintf('$%s$',labs{jnd}),'Interpreter','latex')
      longticks
      xlim([0.9 100])
      set(gca,'XScale','log','XTick',[0 1 10 20 40 100],...
         'YLim',ayl{jnd},'XMinorTick','off')
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
     for ind=1:numel(missarr)
       ftx0=fftshift(fft2(tapersv(:,:,ind)));
       ftph0=imag(ftx0);
       ftph(:,:,ind)=ftph0;
       ftx0=decibel(abs(ftx0).^2);
       % Decibel will produce -Inf if an element is 0; replace these with the
       % next smallest value
       ftx0(isinf(ftx0))=NaN;
       ftx0(isnan(ftx0))=min(ftx0,[],'all');
       ftx(:,:,ind)=ftx0;
       % And again, but this time with padding after conversation w Arthur 4/30/25
       ftx0p=fftshift(fft2(tapersv(:,:,ind),2^nextpow2(params.NyNx(2)),...
                   2^nextpow2(params.NyNx(1))));
       ftph0p=imag(ftx0p);
       ftphp(:,:,ind)=ftph0p;
       ftx0p=decibel(abs(ftx0p).^2);
       % Decibel will produce -Inf if an element is 0; replace these with the
       % next smallest value
       ftx0p(isinf(ftx0p))=NaN;
       ftx0p(isnan(ftx0p))=min(ftx0p,[],'all');
       ftxp(:,:,ind)=ftx0p;
     end
     % Plot with padding
     upad=0;
     if upad; ftp=ftxp; ftphs=ftphp; else; ftp=ftx; ftphs=ftph; end
     for pnd=1:2
     figure()
     numrws=ceil(sqrt(numel(missarr)));
     [ah,ha,H]=krijetem(subnum(numrws,2*ceil(numel(missarr)./numrws)));
     for ind=1:size(ha,2)/2
       % Figures in x-space are in odd columns
       xnd=ind+numrws*floor((ind-1)/numrws);
       % Figures in k-space are in even columns
       knd=ind+numrws*(floor((ind-1)/numrws)+1);
       if ind<=numel(missarr)
         % Display the taper in grid space
         axes(ha(xnd))
         imagesc(tapersv(:,:,ind))
         caxis([0 1])
         hold on
         % Show the decorrelation length scale (the correlation is 1/3 of the
         % variance at a distance of pi*rho) 
         dc(ind,1:2)=circ(pi*th1(3),[],[pi*th1(3) params.NyNx(2)-pi*th1(3)]); 
         set(dc(ind,1),'color','b'); delete(dc(ind,2))
         title(sprintf('$%0.1f\\%%$: Tx;',100*prcmiss(ind)),'Interpreter','latex')
         longticks
         axis image
         shrink(ha(xnd),1.05,1.05)
         set(ha(xnd),'XTick',[1 params.NyNx(2)],'YTick',[1 params.NyNx(1)])
         xlabel('$x_1$','Interpreter','latex')
         ylabel('$x_2$','Interpreter','latex')
         % Display the taper in wavenumber space
         axes(ha(knd))
         colormap(ha(knd),parula)
         if pnd==1
           % square modulus pf the FT of the taper
           imagesc(ftp(:,:,ind))
           caxis(minmax(ftp(:,:,2:end)))
           title('$|FFT2(Tx)|^2$','Interpreter','latex')
         else
           % Phase of the FT of the taper
           imagesc(ftphs(:,:,ind))
           %caxis(minmax(ftphs(:,:,2:end)))
         end

         longticks
         axis image
         shrink(ha(knd),1.05,1.05)
         if upad==1
           set(ha(knd),'XTick',[1 params.NyNx(2)],'YTick',[1 params.NyNx(1)])
         else
           set(ha(knd),'XTick',[1 size(ftxp(:,:,1),1)],'YTick',[1 size(ftxp(:,:,1),2)])
         end
         xlabel('$k_1$','Interpreter','latex')
         ylabel('$k_2$','Interpreter','latex')
         set(ha(knd),'YAxisLocation','right')

         moveh([ha(xnd) ha(knd)],-0.01*(3-4*floor((ind-1)/numrws)))
         moveh(ha(knd),-0.02)

         colormap(ha(xnd),gray(2))
       else
         delete([ha(xnd) ha(knd)])
       end
     end
     cbx=colorbar(ha(numrws),'Location','southoutside');
     cbx.Position(1)=0.15;
     cbx.Position(2)=0.03;
     cbx.Position(3)=0.3;
     cbx.Position(4)=0.02;
     cbx.Label.String='spatial grid space';
     cbx.Label.Position(2)=2.3;
     cbx.XTick=[0 1];
     cbk=colorbar(ha(2*numrws),'Location','southoutside');
     cbk.Position(1)=0.6;
     cbk.Position(2)=0.03;
     cbk.Position(3)=0.3;
     cbk.Position(4)=0.02;
     if pnd==1
       cbk.Label.String='wavenumber space: amplitude (10log10(\cdot))';
     else
       cbk.Label.String='wavenumber space: phase';
     end
     cbk.Label.Position(2)=2.3;

     sti(2)=sgtitle(...
       sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s%i%s%i%s%s',...
       '\sigma^2=',th1(1),' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),...
       ' m; ',params.NyNx(1).*params.dydx(1),' m x ',...
       params.NyNx(2).*params.dydx(2),'m; dy=',params.dydx(1),'m, dx=',...
       params.dydx(2),'m; ',upper(masknm)),...
      'Interpreter','latex');
     % Pause to confirm the figure is nice before saving
     saveas(gcf,...                                                               
       sprintf('covgammiosl_demo8_missing_Ny%iNx%i_dy%idx%i_s2%inu%irh%i_%s_2_%i.eps',...
         params.NyNx,params.dydx,th1.*[1 10 1],date,pnd),'epsc') 
     end
   else
     % flg==1
     % this is the publication style figure where we consider the random case
     % and load three datasets with distinct values of rho
     % params.mask='random';
     % datum='30-Apr-2025';% for random
     % datum='03-Dec-2025';% for random
     datum='20-Jan-2026';
     rho1=2; % rho=2
     rho2=4; % rho=4
     rho3=8; % rho=8
     ds1=load(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',rho1,datum));
     ds2=load(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',rho2,datum));
     ds3=load(sprintf('covgammiosl_demo8vars_flg1_rho%i_%s.mat',rho3,datum));

     % A few things we will need for the plot
     numreals=size(ds1.thhats,1);
     labs2={'$\sigma^2$','$\nu$','$\rho$'};
     labs={labs2{:},'\sigma^2,\nu','\sigma^2,\rho','\nu,\rho',labs2{:}};
     clrs=["#0033FF","#3399FF","#33ccFF"];
     % for consistency between experiments that we are studying, set same limits
     ayl={[0.1 1.2],[0.000001 0.03],[0.4 2.5],[0.001 0.015],[0.1 2.1],[0.005 0.08]};
     covl='$|\mathrm{cov}\{\theta_i,\theta_j\}|$';
     biasl='$(\langle\hat{\theta}\rangle-\theta_0)/\theta_0$';

     for ind=1:np
       ds1.bias(:,ind)=squeeze(nanmean(ds1.thhats(:,ind,:)-ds1.th1(ind)));
       ds1.relbsme(:,ind)=ds1.bias(:,ind)/ds1.th1(ind);
       ds1.relbssd(:,ind)=squeeze(sqrt(var(ds1.thhats(:,ind,:)-ds1.th1(ind),"omitnan"))/...
                       ds1.th1(ind));
       ds1.rmse(:,ind)=sqrt(squeeze(var(ds1.thhats(:,ind,:)))+ds1.bias(:,ind).^2);
       %
       ds2.bias(:,ind)=squeeze(nanmean(ds2.thhats(:,ind,:)-ds2.th1(ind)));
       ds2.relbsme(:,ind)=ds2.bias(:,ind)/ds2.th1(ind);
       ds2.relbssd(:,ind)=squeeze(sqrt(var(ds2.thhats(:,ind,:)-ds2.th1(ind),"omitnan"))/...
                       ds2.th1(ind));
       ds2.rmse(:,ind)=sqrt(squeeze(var(ds2.thhats(:,ind,:)))+ds2.bias(:,ind).^2);
       %
       ds3.bias(:,ind)=squeeze(nanmean(ds3.thhats(:,ind,:)-ds3.th1(ind)));
       ds3.relbsme(:,ind)=ds3.bias(:,ind)/ds3.th1(ind);
       ds3.relbssd(:,ind)=squeeze(sqrt(var(ds3.thhats(:,ind,:)-ds3.th1(ind),"omitnan"))/...
                       ds3.th1(ind));
       ds3.rmse(:,ind)=sqrt(squeeze(var(ds3.thhats(:,ind,:),"omitnan"))+ds3.bias(:,ind).^2);
     end

     % Take care of the rest of the calculations 
     % Prepare a vector for making predictions for the best fit lines in this
     % strange 10log10 v log2 space
     gv=0.8*glens(1):0.1:1.2*glens(end);
     log2sqK=log2(Ks);
     gv=0.98*log2sqK(1):0.1:1.02*log2sqK(end);
     for jnd=1:np
       % first dataset
       mbfrmse1(jnd,:)=polyfit(log2sqK,10*log10(abs(ds1.rmse(:,jnd))),1);
       predrmse1(jnd,:)=polyval(mbfrmse1(jnd,:),gv);
       % second dataset
       mbfrmse2(jnd,:)=polyfit(log2sqK,10*log10(abs(ds2.rmse(:,jnd))),1);
       predrmse2(jnd,:)=polyval(mbfrmse2(jnd,:),gv);
       % third dataset
       mbfrmse3(jnd,:)=polyfit(log2sqK,10*log10(abs(ds3.rmse(:,jnd))),1);
       predrmse3(jnd,:)=polyval(mbfrmse3(jnd,:),gv);
     end
     for jnd=1:np*2
       % first dataset
       mbfemp1(jnd,:)=polyfit(log2sqK(1:end),10*log10(abs(ds1.covet(1:end,jnd))),1);
       mbfcov1(jnd,:)=polyfit(log2sqK(1:end),10*log10(abs(ds1.cov3t(1:end,jnd))),1);
       predemp1(jnd,:)=polyval(mbfemp1(jnd,:),gv);
       predcov1(jnd,:)=polyval(mbfcov1(jnd,:),gv);
       % second dataset
       mbfemp2(jnd,:)=polyfit(log2sqK(1:end),10*log10(abs(ds2.covet(1:end,jnd))),1);
       mbfcov2(jnd,:)=polyfit(log2sqK(1:end),10*log10(abs(ds2.cov3t(1:end,jnd))),1);
       predemp2(jnd,:)=polyval(mbfemp2(jnd,:),gv);
       predcov2(jnd,:)=polyval(mbfcov2(jnd,:),gv);
       % third dataset
       mbfemp3(jnd,:)=polyfit(log2sqK(1:end),10*log10(abs(ds3.covet(1:end,jnd))),1);
       mbfcov3(jnd,:)=polyfit(log2sqK(1:end),10*log10(abs(ds3.cov3t(1:end,jnd))),1);
       predemp3(jnd,:)=polyval(mbfemp3(jnd,:),gv);
       predcov3(jnd,:)=polyval(mbfcov3(jnd,:),gv);
     end

     % Set up the figure
     figure(1);
     clf
     figip=get(gcf,'InnerPosition');
     set(gcf,'InnerPosition',[figip(1:2) 950 1200])
     [ah,ha,H]=krijetem(subnum(6,3));
     % Top row: plot the relative bias of the estimates from the experiments
     % with errorbars displaying the standard deviation. This set of y-axes will
     % be in a linear scale reporting the actual relative bias.
     for ind=1:np
       rnd=ind;
       axes(ah(rnd));
       pl3(ind)=errorbar(log2sqK,ds3.relbsme(:,ind),ds3.relbssd(:,ind),'o',...
                'MarkerEdgeColor',clrs{3},'MarkerFaceColor',clrs{3},...
                'LineWidth',1,'CapSize',5,'Color',clrs{3},...
                'DisplayName',...
                sprintf('$\\rho$=%i: obs.',ds3.th1(3)));
       hold on
       pl2(ind)=errorbar(log2sqK,ds2.relbsme(:,ind),ds2.relbssd(:,ind),'o',...
                'MarkerEdgeColor',clrs{2},'MarkerFaceColor',clrs{2},...
                'LineWidth',1,'CapSize',5,'Color',clrs{2},...
                'DisplayName',...
                sprintf('$\\rho$=%i: obs.',ds2.th1(3)));
       pl1(ind)=errorbar(log2sqK,ds1.relbsme(:,ind),ds1.relbssd(:,ind),'o',...
                'MarkerEdgeColor',clrs{1},'MarkerFaceColor',clrs{1},...
                'LineWidth',1,'CapSize',5,'Color',clrs{1},...
                'DisplayName',...
                sprintf('$\\rho$=%i: obs.',ds1.th1(3)));
       box on
       % Legend
       [leg0(ind),legic]=legend('interpreter','latex','EdgeColor','none',...
                                'box','off','NumColumns',1);
       leg0(ind).AutoUpdate='off';
       % Have the legend markers take up less space
       movelcon(legic)
       % Add the 0-line: estimate equal to the truth
       yline(0)
       if ind==1
         ylabel('$(\langle\hat{\theta}\rangle-\theta_0)/\theta_0$','interpreter','latex');
       end
       ti(rnd)=title(labs2{ind},'Interpreter','latex');
       longticks
       if ind~=1
         set(gca,'YTickLabel','')
       end
       xlim(minmax(log2sqK).*[0.99 1.01])
     end
     for jnd=1:9
       rnd=jnd+np;
       axes(ah(rnd))
       hold on
       if jnd<7
         % Second and third row: the parameter (co)variance with 10log10 scaling
         % Plot the observed values first
         % The first dataset
         plot(log2sqK,10*log10(abs(ds1.covet(:,jnd))),'o',...
              'MarkerEdgeColor',clrs{1},'MarkerFaceColor',clrs{1},...
              'DisplayName',...
              sprintf('$\\rho$=%i: obs. %0.1f',ds1.th1(3),abs(mbfemp1(jnd,1))));
         plot(gv,predemp1(jnd,:),...
                 'Color',clrs{1},'LineStyle','-',...
                 'HandleVisibility','off');
         % The second dataset
         plot(log2sqK,10*log10(abs(ds2.covet(:,jnd))),'o',...
              'MarkerEdgeColor',clrs{2},'MarkerFaceColor',clrs{2},...
              'DisplayName',...
              sprintf('$\\rho$=%i: obs. %0.1f',ds2.th1(3),abs(mbfemp2(jnd,1))));
         plot(gv,predemp2(jnd,:),...
                 'Color',clrs{2},'LineStyle','-',...
                 'HandleVisibility','off');
         % The third dataset
         plot(log2sqK,10*log10(abs(ds3.covet(:,jnd))),'o',...
              'MarkerEdgeColor',clrs{3},'MarkerFaceColor',clrs{3},...
              'DisplayName',...
              sprintf('$\\rho$=%i: obs. %0.1f',ds3.th1(3),abs(mbfemp3(jnd,1))));
         plot(gv,predemp3(jnd,:),...
                 'Color',clrs{3},'LineStyle','-',...
                 'HandleVisibility','off');

         % Plot the predicted values second 
         % The first dataset
         plot(log2sqK,10*log10(abs(ds3.cov3t(:,jnd))),'^',...
              'MarkerEdgeColor',clrs{1},'MarkerFaceColor',clrs{1},...
              'DisplayName',...
              sprintf('pred. %0.1f',abs(mbfcov1(jnd,1))));
         plot(gv,predcov1(jnd,:),...
                'Color',clrs{1},'LineStyle','-','HandleVisibility','off');
         % The second dataset
         plot(log2sqK,10*log10(abs(ds2.cov3t(:,jnd))),'^',...
              'MarkerEdgeColor',clrs{2},'MarkerFaceColor',clrs{2},...
              'DisplayName',...
              sprintf('pred. %0.1f',abs(mbfcov2(jnd,1))));
         plot(gv,predcov2(jnd,:),...
             'Color',clrs{2},'LineStyle','-','HandleVisibility','off');
         % The third dataset
         plot(log2sqK,10*log10(abs(ds3.cov3t(:,jnd))),'^',...
              'MarkerEdgeColor',clrs{3},'MarkerFaceColor',clrs{3},...
              'DisplayName',...
              sprintf('pred. %0.1f',abs(mbfcov3(jnd,1))));
         plot(gv,predcov3(jnd,:),...
              'Color',clrs{3},'LineStyle','-','HandleVisibility','off');
         xlim(minmax(glens).*[0.99 1.01])
       end
       % Fourth row: the parameter rmse with 10log10 scaling
       if jnd>6
         rnd=jnd-6;
         % Plot the observed values first
         % The first dataset
         plot(log2sqK,10*log10(abs(ds1.rmse(:,rnd))),'o',...
              'MarkerEdgeColor',clrs{1},'MarkerFaceColor',clrs{1},...
              'DisplayName',...
              sprintf('$\\rho$=%i: obs. %0.1f',ds1.th1(3),abs(mbfrmse1(rnd,1))));
         plot(gv,predrmse1(rnd,:),...
             'Color',clrs{1},'LineStyle','-','HandleVisibility','off');
         % The second dataset
         plot(log2sqK,10*log10(abs(ds2.rmse(:,rnd))),'o',...
              'MarkerEdgeColor',clrs{2},'MarkerFaceColor',clrs{2},...
              'DisplayName',...
              sprintf('$\\rho$=%i: obs. %0.1f',ds2.th1(3),abs(mbfrmse2(rnd,1))));
         plot(gv,predrmse2(rnd,:),...
              'Color',clrs{2},'LineStyle','-','HandleVisibility','off');
         % The third dataset
         plot(log2sqK,10*log10(abs(ds3.rmse(:,rnd))),'o',...
              'MarkerEdgeColor',clrs{3},'MarkerFaceColor',clrs{3},...
              'DisplayName',...
              sprintf('$\\rho$=%i: obs. %0.1f',ds3.th1(3),abs(mbfrmse3(rnd,1))));
         plot(gv,predrmse3(rnd,:),...
             'Color',clrs{3},'LineStyle','-','HandleVisibility','off');
         xlim(minmax(glens).*[0.99 1.01])
       end
     end

     %%
     % Set consistent axes and plot reference lines for the first row
     allylims=[ylim(ah([1:3]))];
     allylims=minmax([allylims{:}]);
     % Make a little more room for the legends
     allylims(2)=diff(allylims)*1.1+allylims(1);
     for jnd=1:3
       axes(ah(jnd))
       ah(jnd).YLim=allylims;
     end
     % Set consistent axes and plot reference lines for the second and third rows
     allylims=[ylim(ah([4:9]))];
     allylims=minmax([allylims{:}]);
     % Make a little more room for the legends
     allylims(2)=diff(allylims)*1.1+allylims(1);
     for jnd=1:6
       mnd=jnd+np;
       axes(ah(mnd))
       % The reference lines
       rnd=0;
       refl=sqrt(10)./sqrt(2.^(2.*gv));
       for ind=allylims(1)-2:5:allylims(end)+30
         rnd=rnd+1;
         refd(rnd)=plot(gv,ind+10*log10(refl),'-',...
                        'Color',[0.85 0.85 0.85],'HandleVisibility','off');
       end
       ah(mnd).YLim=allylims;
     end
     % Set consistent axes and plot reference lines for the fourth row
     allylims=[ylim(ah([10:12]))];
     allylims=minmax([allylims{:}]);
     % Make a little more room for the legends
     allylims(2)=diff(allylims)*1.9+allylims(1);
     for jnd=10:12
       axes(ah(jnd))
       % The reference lines
       rnd=0;
       refl=sqrt(10)./sqrt(2.^(2.*gv));
       for ind=allylims(1)-2:5:allylims(end)+40
         rnd=rnd+1;
         refd(rnd)=plot(gv,ind+10*log10(refl),'-',...
                        'Color',[0.85 0.85 0.85],'HandleVisibility','off');
       end
       ah(jnd).YLim=allylims;
       xlim(minmax(glens).*[0.99 1.01])
     end
     %%
     % these were the values in the linear space for the (co)variance plots
     % ayl={[0.1 1.2],[0.000001 0.03],[0.4 2.5],[0.001 0.015],[0.1 2.1],[0.005 0.08]};
     %%

     % Set the y-axis labels for rows 2--4
     for jnd=1:9
       rnd=np+jnd;
       axes(ah(rnd))
       if jnd==1
          ylabel(sprintf('$%s%s$',...
          '10\times \mathrm{log}_{10}',...
          '(|\mathrm{var}\{\theta_i\}|'),...
          'Interpreter','latex')
       elseif jnd==4
          ylabel(sprintf('$%s%s$',...
          '10\times \mathrm{log}_{10}',...
          '|\mathrm{cov}\{\theta_i,\theta_j\}|'),...
          'Interpreter','latex')
       elseif jnd==7
          ylabel(sprintf('$%s%s$',...
          '10\times \mathrm{log}_{10}',...
          '|\mathrm{rmse}\{\theta_i\}|'),...
          'Interpreter','latex')
       end
       ti(rnd)=title(sprintf('$%s$',labs{jnd}),'Interpreter','latex');
       longticks
       box on
       if ~ismember(jnd,[1 4 7])
         set(gca,'YTickLabel','')
       end
       if jnd<7
         leg1(jnd)=legend('interpreter','latex','EdgeColor','none',...
                          'box','on','NumColumns',2,'FontSize',8);
         leg1(jnd).AutoUpdate='off';
       else
         leg1(jnd)=legend('interpreter','latex','EdgeColor','none',...
                          'box','on','NumColumns',1,'FontSize',8);
         leg1(jnd).AutoUpdate='off';
       end
     end      
     % Set the x-axis labels
     for jnd=1:12
       xtks=log2sqK;
       xtklbs=string(floor(Ks));
       %reshape(round(sqrt(sum(ds1.tapersv,[1 2]))),[1 6 1]));
       set(ah(jnd),'XTick',xtks,'XTickLabel',xtklbs,'XMinorTick','off','YMinorTick','off')
       if jnd>9
         ah(jnd).XLabel.String='$\sqrt{K}$';
         ah(jnd).XLabel.Interpreter='latex';
       end
       %%
       % Create a second x-axis that quotes the percent of the grid that is 
       % observed for the plotted trial
       xlb='Observed grid (%)';
       xxtks=log2sqK;
       xxtklbs=string(prcobs);
       if jnd<4
         [axx(jnd),xl(jnd),yl(jnd)]=xtraxis(ah(jnd),xtks,xxtklbs,...
           xlb);
       else
         [axx(jnd),xl(jnd),yl(jnd)]=xtraxis(ah(jnd),xtks,xxtklbs,...
           []);
       end
       longticks
       box on
       ah(jnd).XRuler.TickLabelGapOffset=-2;
       ah(jnd).XTickLabelRotation=0;
       axx(jnd).XRuler.TickLabelGapOffset=-3;
    end
  
    byl=[-0.8 0.8];
    for ind=1:np
      axes(ah(ind))
      set(gca,'YLim',byl)
      box on
    end
    for ind=1:6
      xnd=12+ind;
      % Display the taper in grid space
      axes(ah(xnd))
      imagesc(ds1.tapersv(:,:,ind))
      caxis([0 1])
      hold on
      % Show the decorrelation length scale (the correlation is 1/3 of the
      % variance at a distance of pi*rho)
      dc3(ind,1:2)=circ(pi*ds3.th1(3),[],...
                        [pi*ds3.th1(3) size(ds3.tapersv,2)-pi*ds3.th1(3)]);
      set(dc3(ind,1),'color',clrs{3},'LineWidth',2); delete(dc3(ind,2))
      dc2(ind,1:2)=circ(pi*ds2.th1(3),[],...
                       [pi*ds2.th1(3) size(ds2.tapersv,2)-pi*ds2.th1(3)]);
      set(dc2(ind,1),'color',clrs{2},'LineWidth',2); delete(dc2(ind,2))
      dc1(ind,1:2)=circ(pi*ds1.th1(3),[],...
                       [pi*ds1.th1(3) size(ds1.tapersv,2)-pi*ds1.th1(3)]);
      set(dc1(ind,1),'color',clrs{1},'LineWidth',2); delete(dc1(ind,2))

      title(sprintf('$%0.1f\\%%$ observed',prcobs(ind)),'Interpreter','latex')
      longticks
      axis image
      shrink(ah(xnd),1.05,1.05)
     % set(ah(xnd),'XTick',[1 params.NyNx(2)])
      set(ah(xnd),'XTick','','YTick','')
      xlabel('$x_1$','Interpreter','latex')
      if ind==1
     %   set(ah(xnd),'YTick',[1 params.NyNx(1)])
        ylabel('$x_2$','Interpreter','latex')
      end

      colormap(ah(xnd),gray(2))
   end
   tprnd=12;
   ahhdiff=(ah(tprnd+1).OuterPosition(1)-ah(tprnd-2).OuterPosition(1))*1.5;
   ah(tprnd+1).OuterPosition(1)=ah(tprnd+1).OuterPosition(1)-ahhdiff;
   for xnd=tprnd+2:tprnd+6
     ah(xnd).OuterPosition(3)=ah(tprnd+1).OuterPosition(3);
     ah(xnd).OuterPosition(4)=ah(tprnd+1).OuterPosition(4);
   end
   ah(tprnd+2).OuterPosition(1)=0.6*ah(tprnd+1).OuterPosition(3)+ah(tprnd+1).OuterPosition(1);
   ah(tprnd+3).OuterPosition(1)=0.6*ah(tprnd+2).OuterPosition(3)+ah(tprnd+2).OuterPosition(1);
   ah(tprnd+4).OuterPosition(1)=0.48*ah(tprnd+3).OuterPosition(3)+ah(tprnd+3).OuterPosition(1);
   ah(tprnd+5).OuterPosition(1)=0.6*ah(tprnd+4).OuterPosition(3)+ah(tprnd+4).OuterPosition(1);
   ah(tprnd+6).OuterPosition(1)=0.6*ah(tprnd+5).OuterPosition(3)+ah(tprnd+5).OuterPosition(1);
   ah(tprnd+4).OuterPosition(2)=ah(tprnd+1).OuterPosition(2);
   ah(tprnd+5).OuterPosition(2)=ah(tprnd+1).OuterPosition(2);
   ah(tprnd+6).OuterPosition(2)=ah(tprnd+1).OuterPosition(2);
   movev(ah(tprnd+1:tprnd+6),-0.02)
   cbx=colorbar();
   cbx.XTick=[0.25 0.75];
   cbx.XTickLabel={'0' '1'};
   moveh(ah(tprnd+1:tprnd+6),-0.05)
   ah(tprnd+6).OuterPosition(3)=ah(tprnd+5).OuterPosition(3);
   ah(tprnd+6).OuterPosition(1)=2*ah(tprnd+5).OuterPosition(1)-ah(tprnd+4).OuterPosition(1);
   moveh(ah([13:18]),0.02)
  
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
   movev([ah axx],-0.03)

   % Some more moving and shrinking
   shrink([ah([1:tprnd]) axx([1:tprnd])],0.9,0.9)
   movev(ah([tprnd+1:tprnd+6]),-0.13)
   movev([ah([10:12]) axx([10:12])],-0.11)
   movev([ah([7:9]) axx([7:9])],-0.08)
   movev([ah([4:6]) axx([4:6])],-0.05)
   movev([ah([1:3]) axx([1:3])],-0.02)
  % moveh(ah([tprnd+4:tprnd+6]),-0.025)
   movev(xl,0.1)
   % move the titles 
   set(ti([1:12]),'Units','pixels')
   movev(ti([1:12]),15)
   for ind=1:3 
     axx(ind).XLabel.Units='pixel';
     axx(ind).XLabel.Position(2)=axx(ind).XLabel.Position(2)+15;
   end
   leg0(1).IconColumnWidth=leg1(7).IconColumnWidth;
   leg0(2).IconColumnWidth=leg1(7).IconColumnWidth;
   leg0(3).IconColumnWidth=leg1(7).IconColumnWidth;
   moveh(leg0,0.01)

   % Final tinkering
   keyboard
   for ind=1:12
     if ind<10
       ah(ind).XTickLabel='';
     end
     axes(axx(ind))
     if ind>3
       axx(ind).XTickLabel='';
     end
   end
   movev(leg0([1:3]),0.006)
   movev([ah([4:12]) axx([4:12])],0.023)
   movev(leg1([1:9]),0.014)
   movev([ah([7:12]) axx([7:12])],0.023)
   movev(leg1([4:9]),0.023)
   movev([ah([10:12]) axx([10:12])],0.023)
   movev(leg1([7:9]),0.023)
   for ind=1:12;ti(ind).FontSize=11;end
   for ind=4:9;ah(ind).YLim=[-40 40];end
   movev(ah([end-5:end]),0.09)

   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo8_missing_Ny%iNx%i_dy%idx%i_%s_final.eps',...
       params.NyNx,params.dydx,date),'epsc')
  end
elseif strcmp(th,'demo9')
   % Growing domain asymptotics study looking at whether the blurred Fisher 
   % inverse evaluated at the truth converges to the parameter covariance 
   % calculated at the truth with increasing sample size; all more idealized
   % than one may expect in practice. Compare to DEMO6, which does consider
   % sampling.
   params=[];params.dydx=[1 1];params.blurs=-1;params.taper=1;ifinv=[1 1 1];
   np=3;

   % Labels
   labs={'\sigma^2', '\nu', '\rho', '\sigma^2,\nu',...
         '\sigma^2,\rho', '\nu,\rho'};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];
   % Font size for latex interpreted text
   fs=11;

   % For now, small values so that we have a chance at forming the data vector
   % and finding an estimate for small grids 
   rh=4;
   th1=[3 0.75 rh];
   % ARTHURS EXAMPLES
   rh=8;
   th1=[1 0.5 rh];
   % Sample growing grid length as multiples of pi rho from logspace 
   glens=2:8;
   % TEMPORARY FOR ARTHURS
   glens=[6 7 8 9 10 10.5];
   log2glens=2.^(glens);

   try
     % If we have already calculated the estimates, load them
     datum='20-Jan-2026';
     % TEMPORARY for Arthur's experiment
    % [1 0.5 8]; 
    datum='20-Mar-2026';
    % [1 10 8]; datum='19-Mar-2026';
    % [1 0.33 8]; datum='18-Mar-2026';
     load(sprintf('covgammiosl_demo9_calcs_%s.mat',datum))
     % The calculated covariance will be exact; calculated using method 2 when
     % the grid size is less than about 120^2, then using method 3 when the grid
     % size is larger. These methods are always the same.
   catch
     for ind=1:size(glens,2)
       % Used to do this
       % Create a larger grid with length linearly increasing in pi rho
       % gs=rh*pi*loggs(ind); 
       % Now the grid size is in terms of Ny=Nx evenly spaced in log2
       gs=log2glens(ind); 
       params.NyNx=[floor(gs) floor(gs)];
       % If we want to look at irregular fields,
       if isfield(params,'mask'); if strcmp(params.mask,'random'); 
           [~,I]=muckit(randn(params.NyNx),params,0.667); 
           tapersv{ind}=I;
           params.taper=tapersv{ind};
       end; end
       % Physical grid size, prod(NyNx / dydx)
       Ns(ind)=prod(params.NyNx./params.dydx);

     % TEMPORARY for Arthur's examples
       covthi=covthosl(th1,params,1,ifinv);

     %  if prod(params.NyNx)<100^2
     %    % Calculate the covariance using the 'dftmtx' method
     %    covthi=covthosl(th1,params,2,ifinv);
     %  else
     %    % Calculate the covariance using the 'diagonals' method
     %    meth=3;

     %    disp(sprintf('Calculating COVTHOSL for meth %i and NyNx %i %i',...
     %            meth,params.NyNx))
     %    covthi=covthosl(th1,params,meth,ifinv);
     %  end
       k=knums(params);
       [~,covFi]=fishiosl(k,th1,params);
       % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
       covth(ind,:)=trilos(covthi);
       covF(ind,:) =trilos(covFi);
     end
     eff=covF./covth;
     % Consider saving your calculations by setting 'svcalcs'
     svcalcs=1;
     if svcalcs
      save(sprintf('covgammiosl_demo9_calcs_%s.mat',datum),...
                   'covth','covF','params','th1','Ns','eff')
     end
   end

   figure(1);
   clf
   figip=get(gcf,'InnerPosition');
   set(gcf,'InnerPosition',[figip(1:2) 950 475])
   [ah,ha,H]=krijetem(subnum(2,3));
   colororder({'k',clrs{3}})
   gv=0.8*glens(1):0.1:1.2*glens(end);
   % The trend is not very linear if we look at sample sizes smaller than this
   gv5=0.8*max(5,glens(1)):0.1:1.2*glens(end);
   for jnd=1:6
     axes(ah(jnd))
     yyaxis left
     hold on
     % Fit a slope to the covariance calculations as a function of grid length
     % in log-log
     mbfcovth(jnd,:) =polyfit(glens(glens>=5),10*log10(abs(covth(glens>=5,jnd))),1);
     mbfcovF(jnd,:)  =polyfit(glens(glens>=5), 10*log10(abs(covF(glens>=5,jnd))),1);
     mbfratio(jnd,:) =polyfit(glens(glens>=5),  10*log10(abs(eff(glens>=5,jnd))),1);
     predcovth(jnd,:)=polyval(mbfcovth(jnd,:),gv5);
     predcovF(jnd,:) =polyval(mbfcovF(jnd,:), gv5);
     predratio(jnd,:)=polyval(mbfratio(jnd,:),gv5);

     % Plot the absolute covariance calculated from covthosl and the fisher
     % inverse
     ll1=plot(glens,10*log10(abs(covth(:,jnd))),...
            'o','MarkerEdgeColor','k',...
            'MarkerFaceColor','k');
     ll2=plot(glens,10*log10(abs(covF(:,jnd))),...
            's','MarkerEdgeColor','k',...
            'MarkerFaceColor','k');
     yyaxis right
     ll3=plot(glens,10*log10(eff(:,jnd)),...
            '*','MarkerEdgeColor',clrs{3},...
            'MarkerFaceColor',clrs{3});
     % Plot the best fit lines
     yyaxis left
     g(1,jnd)=plot(gv5,predcovth(jnd,:),'k-','HandleVisibility','off');
     g(2,jnd)=plot(gv5,predcovF(jnd,:),'k-','HandleVisibility','off');
     yyaxis right
     hold on
     g(3,jnd)=plot(gv5,predratio(jnd,:),'Color',clrs{3},'LineStyle','-',...
             'HandleVisibility','off');
   end
   for jnd=1:6
     axes(ah(jnd));
     yyaxis left
     % List slopes in the legend and explain in a figure caption
     [leg1(jnd),legic]=legend(...
       {append('$\mathrm{cov}\{\theta_0\}, ',num2str(round((mbfcovth(jnd,1)),1)),'$'),...
        append('$\bar{\mathcal{F}}^{-1}(\theta_0), ',num2str(round((mbfcovF(jnd,1)),1)),'$'),...
        append('$e(\theta_0)$, $',num2str(round(mbfratio(jnd,1),1)),'$')},...
        'Interpreter','latex','box','off',...
        'FontSize',fs);
     % Have the legend markers take up less space
     movelcon(legic,[0.25 0.15 0.35])
     % Labels, titles, ticks, common limits
     if jnd==1
        ylabel('$10\times \mathrm{log}_{10}(|\mathrm{var}\{\theta_i\}|)$',...
               'Interpreter','latex','FontSize',fs)
     elseif jnd==4
        ylabel('$10\times \mathrm{log}_{10}(|\mathrm{cov}\{\theta_i,\theta_j\}|$)',...
               'Interpreter','latex','FontSize',fs)
     end
     t(jnd)=title(sprintf('$%s$',labs{jnd}),'Interpreter','latex','FontSize',fs);
     longticks
     ah(jnd).XRuler.TickLabelGapOffset=-2;
     box on
     xlim(minmax(glens).*[0.99 1.01])
     xtks=glens;
     xtklbs=string(log2glens);
     set(ah(jnd),...
         'XTick',xtks,'XTickLabel',xtklbs,'XMinorTick','off',...
         'YMinorTick','off')
     xla='$Ny \times Nx$ Samples';
     xla='$\sqrt{K}$';
     if jnd>3
       xlabel(xla,...
              'Interpreter','latex','FontSize',fs)
     end
   end
   % Set consistent axes
   allylims=[ylim(ah([1:6]))];
   allylims=minmax([allylims{:}]);
   % Make a little more room for the legends
   allylims(2)=diff(allylims)*1.5+allylims(1);
   for jnd=1:6
     axes(ah(jnd))
     % Plot reference diagonals as light gray lines on bottom
     rnd=0;
     refl=sqrt(10)./sqrt(2.^(2*(gv)));
     for ind=allylims(1)-2:5:allylims(end)+15
       rnd=rnd+1;
       refd(rnd)=plot(gv,ind+10*log10(refl),'-',...
            'Color',[0.85 0.85 0.85],'HandleVisibility','off');
     end
     ah(jnd).YLim=allylims;
   end

   for jnd=1:6
     axes(ah(jnd))
     yyaxis right
     set(gca,'YMinorTick','off')
     if jnd==3|jnd==6
       ylabel('$10\times \mathrm{log}_{10}(e (\theta_0))$',...
              'Interpreter','latex','FontSize',fs)
     end
     % Create a second x-axis that quotes the length of one side of the grid
     % in units of pi rho
     xlb='Length [$\pi\rho$]';
     xxv=round(linspace(1,ceil(log2glens(end)./(pi*th1(3))),numel(glens)));
     xxv=ceil(log2glens./(pi*th1(3)));
     xxtks=log2(pi*th1(3).*xxv);
     xxtklbs=string(xxv);
     if jnd<4
       [axx(jnd),xl(jnd),yl(jnd)]=...
         xtraxis(ah(jnd),xxtks,xxtklbs,xlb);
       xl(jnd).Interpreter='latex';
       xl(jnd).Units='pixels';
     else
       [axx(jnd),xl(jnd),yl(jnd)]=...
         xtraxis(ah(jnd),xxtks,xxtklbs,[]);
     end
     longticks
     ah(jnd).XRuler.TickLabelGapOffset=-2;
     ah(jnd).XTickLabelRotation=0;
     axx(jnd).XRuler.TickLabelGapOffset=-3;
   end      
   % Set consistent axes
   allylims=[ylim(ah([1:6]))];
   allylims=minmax([allylims{:}]);
   % Make a little more room for the legends
   allylims(2)=diff(allylims)*1.5+allylims(1);
   for jnd=1:6
     ah(jnd).YLim=allylims;
   end
   movev(t,7)
   movev(xl,20)

   % Super title to display the parameters
   ti=sgtitle(sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%0.1f%s%0.1f',...
                  '\sigma^2=',th1(1),' \mathrm{[unit]}^2, \nu=',th1(2),...
                  ', \rho=',th1(3),' m $|$ dy=',params.dydx(1),...
                  ', dx=',params.dydx(2)),...
              'Interpreter','latex');
   % Shift the figure objects to make space for the title
   movev([ah(1:3) axx(1:3)],-0.07)
   movev([ah(4:6) axx(4:6)],-0.01)
   movev([leg1([1:6])],0.005)

   keyboard
   % Final tinkering
   for ind=1:6
     if ind<4
       ah(ind).XTickLabel='';
     end
     axes(axx(ind))
     if ind>3
       axx(ind).XTickLabel='';
     end
   end
   movev([ah([4:6]) axx([4:6]) leg1([4:6])],-0.01)
   for ind=1:6;
     t(ind).FontSize=11;
     t(ind).Units="normalized";
     xl(ind).Units="normalized";
   end
   movev(xl([1:3]),-0.04)
   movev(t([1:3]),0.04)
   delete(ti)
  % ti.FontSize=9;
  % movev([ah([1:end]) axx([1:end]) leg1([1:end])],-0.03)

   % Save the figure
   saveas(gcf,...
     sprintf('covgammiosl_demo9_gridsizefish_s2%inu%irh%i_%s.eps',...
       th1.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo10')
  % What is the cutoff in the smoothness parameter for which the inverse fisher
  % approaches the calculated parameter variance. We must know.
  
  % a large enough grid:
  gs=[1024 floor(2^10.5)];;
  th=[1 0.25 8]; % arthur's s2 and rh
  nus=0.25:0.25:10;
  datum=date;

  try
    load(sprintf('covgammiosl_demo10_calcs_%s.mat',datum))
  catch
    for jnd=1:numel(gs)
      p=[];p.dydx=[1 1];p.NyNx=[gs(jnd) gs(jnd)];p.blurs=-1;
      for ind=1:numel(nus)
        th1(ind,:)=th;th1(ind,2)=nus(ind);
        covthi=covthosl(th1(ind,:),p,1,[1 1 1]);
        k=knums(p);
        [~,covFi]=fishiosl(k,th1(ind,:),p);
        % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
        covth{jnd}(ind,:)=trilos(covthi);
        covF{jnd}(ind,:) =trilos(covFi);
      end
      eff{jnd}=covF{jnd}./covth{jnd};
    end
    % Consider saving your calculations by setting 'svcalcs'
    svcalcs=1;
    if svcalcs
      keyboard
      save(sprintf('covgammiosl_demo10_calcs_%s.mat',datum),...
                   'covth','covF','p','th1','eff')
    end
  end
  figure()
  plot(nus,log10(eff{1}(:,1))-log10(eff{2}(:,1)),'bo')
  hold on
  plot(nus,log10(eff{1}(:,2))-log10(eff{2}(:,2)),'k^')
  plot(nus,log10(eff{1}(:,3))-log10(eff{2}(:,3)),'r*')
  yline(0)
  lg=legend('\sigma^2','\nu','\rho');
  lg.Location='northwest';
  xlabel('\nu_0')
  ylabel('log10(e_{K=1024}) - log10(e_{K=1448})')
  xline(3)
  set(gcf().Children(2).Children(1),'HandleVisibility','off')
  keyboard
  saveas(gcf, sprintf('covgammiosl_demo10_%s.eps',date))
end
