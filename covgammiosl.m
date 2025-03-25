function varargout=covgammiosl(th,params,method,ifinv);
% [covg,grads]=COVGAMMIOSL(th,params,method,ifinv);
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
% th        The Matern parameters, [s2 nu rh].
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
% grads     Specific only to method 1: the score calculated for each of the
%           sampled empirical periodograms; used for demos/validation only.
%
% SEE ALSO:
%
% COVTHOSL, MATERNOSY, BLUROSY
%
% EXAMPLES:
% 
% th=[1 0.5 2]; params=[]; params.NyNx=[8 8]; params.dydx=[1 1];
% params.blurs=Inf; params.taper=1; ifinv=[1 1 1];
% tic; [covg1,grads]=covgammiosl(th,params,1,ifinv); toc
% tic;  covg2       =covgammiosl(th,params,2,ifinv); toc
% tic;  covg3       =covgammiosl(th,params,3,ifinv); toc
% difer(covg2-covg3,5)
% difer(covg1-covg3,1)
%
% % For studying random deletions and other non-unit tapers
%
% th=[1 0.5 1];params=[];params.NyNx=[3 3];params.dydx=[1 1];
% params.blurs=Inf;ifinv=[1 1 1];tpr=ones(params.NyNx);tpr(2,3)=0;
% params.taper=tpr;
% covg2 = covgammiosl(th,params,2,ifinv);
% covg3 = covgammiosl(th,params,3,ifinv);
%
% Last modified by fjsimons-at-princeton.edu, 03/25/2025
% Last modified by olwalbert-at-princeton.edu, 03/25/2025

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

  % Initialize covg
  covg=zeros(np,np);

  % Calculating the score requires the expected blurred periodogram, which we
  % form from BLUROSY through MATERNOSP, which we must not provide blurs=Inf to
  paramsSbar=params;
  paramsSbar.blurs=-1;
  Sbar=maternosp(th,paramsSbar);

  % We also need blurred partial derivatives of the expected periodogram, which
  % are calculated in BLUROSY by calling MATERNOSY; initialize [prod(NyNx) np]
  dSbardth=zeros([size(Sbar) np]);
  for thnd=1:np
    if ifinv(thnd)
        % Only calculate partials for the parameters requested by ifinv
        dSbardthnd=blurosy(th,params,xver,'efs',[],thnd);
    end
    % Store as vectors
    dSbardth(:,:,thnd)=dSbardthnd;
  end

  % Calculate the covariance of the score for the selected method
  switch method
    case 1
      % Calculate the covariance of the score by SAMPLING many realizations of
      % the empirical periodogram

      % Set the number of realizations. 1e3-2e3 is sufficient for terms to
      % asymptote out; if not runinpar, 1e3 takes 8.4s for a 128 by 128 grid, 
      % if runinpar with 23 workers, 1e3 takes 1.1s, with more scaling linearly
      numreals=1000; 
      % Pre-allocate space for each score
      grads=zeros(numreals,3);

      % Simulate repetitively for the same parameters, returning the empirical 
      % periodogram each time
      runinpar=1;
      if runinpar==1
        % Create the parallel pool with NumWorkers; running into problems using
        % all physical cores, despite twice as many logical cores detected
        NumWorkers=feature('numcores')-1;
        if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end
        parfor ind=1:numreals
          % Generate an empirical (modified, tapered, windowed) periodogram
          % through SGP simulation; SIMULOSL already takes the taper in to
          % account and normalizes Semp appropriately
          [~,~,~,~,~,Semp]=simulosl(th,params,xver);
          Semps(:,ind)=Semp;
        end
      else
        Semps=zeros(numel(Sbar),numreals);
        for ind=1:numreals
          [~,~,~,~,~,Semp]=simulosl(th,params,xver);
          Semps(:,ind)=Semp;
        end
      end

      % Calculate the elements of the score for each parameter vectorized in
      % dSbardth [numreals by np] % remember dSbardth./Sbar.^2 <- MAOSL/Sbar
      % (A.53) with Sbar factored out; (A.54) would have gone via cov(semps)
      % grads=deal(squeeze(sum((1-Semps/Sbar).*mths)));
      % grads=deal(squeeze(sum((1-Semps/Sbar).*dSbardth./Sbar)));
      grads=deal(squeeze(sum((Sbar-Semps).*dSbardth./Sbar.^2))); 

      if isfield(params,'taper')&numel(params.taper)==prod(params.NyNx)
         % Adjust for the size of an explicit taper as in SIMULOSL and MLEOSL 
         % Accounting for dydx seems... off.
         normfact=(prod(params.NyNx)./sum(Tx(:).^2)).^2.*prod(params.NyNx).^2;
      else
         normfact=prod(params.NyNx).^2;
      end
      % %% OLW - dy~=dx leads to an inconsistency with methods 2 and 3; 
      % I think that dydx is being applied consistently to Sbar, Semp, and 
      % dSbardth, and it doesn't make sense to scale the score by dydx. Are both
      % methods 2 and 3 off? Should dydx be incorporated in spatmat?
      % Straight calculation of the covariance of the score
      covg=cov(grads)./normfact;
    case 2
      % Exact calculation of the gradient of the empirical periodograms using
      % Isserlis Rule, which allows us to instead calculate the sum of two 2D-FFTs
      % (one as a Hermitian conjugate) of the Matern autocovariance calculated 
      % for a NyNx^2 by NyNx^2 distance grid
      % The largest grid size I can calculate using R2024a on Aguilonius is
      % about 128x128; killed at 140x140

      % Isserlis' Rule allows us to calculate the covariance of the empirical
      % periodogram as a combination of transposed and conjugated Fourier 
      % Transforms of the autocovariance (i.e., the expectation of the data
      % vector; see Walden+1999 eqs 7-8):
      % cov{|Hk|^2,|Hk'|^2} = U'*Cmn*U + U'*Cmn*U.', for the 1-D case

      % The distance grid, dxxp, of all pairwise L2 distances between two common 
      % gridded periodograms that takes into account the spacing in x and y
      dydx=params.dydx;NyNx=params.NyNx;
      ys=0:NyNx(1)-1;xs=0:NyNx(2)-1;
      [grdx,grdy]=meshgrid(xs(:),ys(:));
      grdx=grdx(:).*dydx(2);
      grdy=grdy(:).*dydx(1);
      dxxp=xxpdist([grdx grdy]);
      % The following is equivalent; is it faster to calculate here?
      % lagx=grdx-grdx';
      % lagy=grdy-grdy';
      % dxxp=sqrt(lagx.^2+lagy.^2);

      % The autocovariance at the lag-distances
      Cmn=maternosy(dxxp,th);
      % %% OLW testing is still in progress.
      % We must apply the taper here; note that previously this method was
      % consistent with method 2 for fully sampled grids without incorporating
      % the autocorrelation sequence in the commented out lines below

     % % Produce the autocorrelation sequence for the typical NyNx grid (eq. 12)
     % for ind=grdy'+1
     %   for jnd=grdx'+1
     %     t(ind,jnd)=sum(sum(Tx(1:NyNx(1)-ind+1,1:NyNx(2)-jnd+1).*...
     %                       (conj(Tx(ind:end,jnd:end)))));
     %   end
     % end
     % t=t.*Tx./sum(sum(Tx.^2));

     % % Flatten as a row vector for construction of a nested-block Toeplitz
     % % (circulant) matrix that corresponds to the distance matrix, dxxp 
     % t=t(:)';
     % for ind=1:NyNx(1)
     %   begi=1+(ind-1)*NyNx(1);
     %   endi=NyNx(1)+(ind-1)*NyNx(1);
     %   tb=toeplitz(t(begi:endi));
     %   nxnxcircblock{ind}=tb;%./sum(sum(tb.^2));
     % end
     % tt = cell2mat(nxnxcircblock(toeplitz(1:NyNx(1))));

     % if numel(params.taper)>1
     %   Cmn=Cmn.*tt;
     % end

     % Seems to be missing a normalization? Should we multiply the
     % cross-covariance by the mask?... Probably not because it is isotropic
     % Cmn=Cmn.*reshape(Tx',1,numel(Tx)).*reshape(Tx',1,numel(Tx)).';%.*...
     %      (sum(Tx(:))./prod(params.NyNx)).^2;

     % Notice that Cmn is symmetric and positive-definite; can we get any
     % improvements by decomposing?
     % chCmn=chol(Cmn);

      % % The DFT matrix (mxm)
      Um=dftmtx(size(Cmn,1)); 
      % % The DFT matrix (nxn)
      Un=dftmtx(size(Cmn,2)); 

      % More notes: 
      % In 1-D, for a vector m x 1
      %   fft(Cm)  =(Cm.'*Um).' + O(eps)
      % In 2-D, for an m x n matrix
      %   fft2(Cmn)=Um*Cmn*Un.' + O(eps)
      %            =tensorprod(Un,tensorprod(Um,Cmn.',2,2).',2,1).'+O(eps)
      %            =tensorprod(Un,tensorprod(Um,Cmn  ,2,1)  ,2,2).'+O(eps)
      % The first term of the covariance of the periodogram via Isselis' rule 
      % (eq. 9 of VoWE; also see Walden+1994)
      % In 1-D, EJJht=Um'*Cm*Um./n;
      % In 2-D, EJJht=(Um*(Um*Cmn*Un.')'*Un.')'/prod(NyNx);
      % We have m*n sets of perturbed autocovariance matrices; we will perform
      % the above operations on each

      % The distance matrix is comprised of Nx [Ny by Ny] blocks in the rows and
      % columns. Reshape Cmn such that there are NyNx ``column vectors'' that
      % are Nx sets of the [Ny by Ny] blocks; dimensions should be Nx by Ny by
      % NyNx
      Cm_n_mn=reshape(Cmn,[NyNx(1) NyNx(2) prod(NyNx)]);
      Um=dftmtx(size(Cm_n_mn,1));
      Un=dftmtx(size(Cm_n_mn,2));

      % This is like taking fft2(Cm_n_mn(:,:,i),[],2) for the ith slice
      EJJht_in =pagetranspose(tensorprod(Un,tensorprod(Um,Cm_n_mn,2,1),2,2));
      % Reshape and reorder dimensions 
      EJJht_in =reshape(reshape(...
                    permute(EJJht_in,[3 2 1]),[prod(NyNx) prod(NyNx)]),...
                    [NyNx(1) NyNx(2) prod(NyNx)]);
      EJJht_out=pagetranspose(...
                conj(tensorprod(Un,tensorprod(Um,conj(EJJht_in),2,1),2,2)));
      % Reshape and reorder dimensions 
      EJJht_out=reshape(permute(permute(EJJht_out,[2 1 3]),[2 3 1]),...
                [prod(NyNx) prod(NyNx)]).';

      % The second term of the covariance of the periodogram via Isselis' rule 
      % (eq. 9 of VoWE; also see Walden+1994)
      % In 1-D, EJJt=Um'*Cm*Um.'./NyNx;
      % In 2-D, EJJht=(Um*(Um*Cmn*Un.')*Un)'/prod(NyNx)
      EJJt_out=pagetranspose(tensorprod(Un,tensorprod(Um,EJJht_in,2,1),2,2));
      % Reshape and reorder dimensions 
      EJJt_out=reshape(permute(permute(EJJt_out,[2 1 3]),[2 3 1]),...
        [prod(NyNx) prod(NyNx)]).';
      
      % This should go into a demo
      % While the first and second terms below sum to similar quantities, we can
      % quickly see wave-number dependent distinctions:
      devplt=0;
      if devplt
        clf 
        ah(1)=subplot(321); imagesc(dxxp); cb1=colorbar; 
        cb1.Label.String='euclidean distance [m]'; 
        title('Distance matrix')
        ylabel('x2-lag [m]'); 
        ah(2)=subplot(322); imagesc(Cmn); cb2=colorbar; 
        cb2.Label.String='covariance [unit^2]'; 
        title('C(||x-x''||)')
        ah(3)=subplot(323);
        imagesc(log(abs(EJJht_out./prod(NyNx)).^2));
        cb3=colorbar;cb3.Label.String='log([|unit|^2/m])';
        title('|(U(UCU'')*U)*|^2')
        ylabel('x2-lag [m]'); 
        ah(4)=subplot(324);
        imagesc(log(abs(EJJt_out./prod(NyNx)).^2));
        cb4=colorbar;cb4.Label.String='log([|unit|^2/m])';
        title('|U(UCU'')U|^2')
        % Ratio of Isserlis' terms for periodogram covariance
        ah(5)=subplot(325);
        imagesc(abs(EJJht_out./prod(NyNx)).^2./abs(EJJt_out./prod(NyNx)).^2);
        cb5=colorbar;
        cb5.Label.String='amplitude';
        title('|(U(UCU'')*U)*|^2 / |U(UCU'')U|^2')
        xlabel('x1-lag [m]'); ylabel('x2-lag [m]'); 
        ah(6)=subplot(326);
        imagesc(abs(EJJt_out./prod(NyNx)).^2./abs(EJJht_out./prod(NyNx)).^2);
        cb6=colorbar;
        cb6.Label.String='amplitude';
        title('|U(UCU'')U|^2/|(U(UCU'')*U)*|^2')
        xlabel('x1-lag [m]');
        for ind=1:size(ah,2)
          ah(ind).XTick=1:numel(grdx);ah(ind).XTickLabel=grdx;
          ah(ind).YTick=1:numel(grdy);ah(ind).YTickLabel=grdy;
          ah(ind).XTickLabelRotation=0;ah(ind).YTickLabelRotation=90;
          axes(ah(ind));longticks; axis image; box on
        end
        ti=sgtitle(sprintf('th=[%0.2g %g %0.2g] | %ix%i',th,params.NyNx));
        movev(ah,-0.02)
        figdisp('covgammiosl',sprintf('%i_%i_%i_%i_%i',round(th),params.NyNx),[],1);
      end

      if isfield(params,'taper')&numel(params.taper)==prod(NyNx)
        normfact=(prod(params.NyNx)./sum(Tx(:).^2)).^2.*prod(params.NyNx).^2;
      else
        normfact=prod(NyNx).^2;
      end

      % Sum and normalize the quantities above; eq. 8 of Walden+1994
      Hk2cov=(abs(EJJht_out).^2+abs(EJJt_out).^2)./prod(NyNx).^2;

      % Un-do the normalization factor applied in BLUROSY for the blurred
      % expected periodogram and its partial derivatives; shift both such
      % that the zero-wavenumber is in the 1,1 position 
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

      % Product of all partials in row-order: 
      % [s2s2 s2nu s2rh nus2 nunu nurh rhs2 rhnu rhrh]
      % f_seq=[f(:,1).*f f(:,2).*f f(:,3).*f];
      % Arthur special ?
      normalize=0;
      if normalize
        Sbar_n=v2s(Sbar,params);
        seq_ep=Sbar_n.*Sbar_n;
      else
        seq_ep=1;
      end 

      % Final assembly (A54)
      covg=fax'*Hk2cov*fax./seq_ep./normfact;
      % I know that I need to multiply by prod(dydx).^2, but I am not sure where
      % would be best. Doing so here works for now.
      covg=covg.*prod(dydx).^2;
      % wkspc=whos;sum([wkspc.bytes])
      grads=NaN;
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
      f=dSbardth./Sbar.^2;

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
            fftcov_diag=blurosy(th,paramsSbar,xver,'efd',[my mx]);                        
            fftcov_diag=v2s(fftcov_diag,params)*nf;                                           
            cdd=fftcov_diag(max([1,my+1]):min([Ny+my,Ny]),...                         
                            max([1,mx+1]):min([Nx+mx,Nx]));                           
            cdd_l2=abs(cdd).^2;                                                  
            cdd_l2=cdd_l2(:);                                                       
                                                                                 
            % Use the same framework to select the antidiagonals from BLUROSY
            fftcov_anti=blurosy(th,paramsSbar,xver,'efd',[my2 mx2]);
            fftcov_anti=v2s(fftcov_anti,params)*nf;
            cad=fftcov_anti(max([1,my2-Ny+2]):min([my2+1,Ny]),...
                            max([1,mx2-Nx+2]):min([mx2+1,Nx]));
            cad_l2=abs(cad).^2;
            cad_l2=cad_l2(:);

            % We need to be a bit more verbose about this in parfor
            if normalize
              seq_ep=Sbar_n(ep_ind1).*Sbar_n(ep_ind2);
              Sbar_n3=Sbar_n(ep_ind3);
              seq_ep2=Sbar_n3.*flip(Sbar_n3,1);
            else
              seq_ep =1;
              seq_ep2=1;
            end

            % The partials indexed for the current `diagonal'
            f1=f(ep_ind1,:);f2=f(ep_ind2,:);
            % Product of all partials in row-order:
            % [s2s2 s2nu s2rh nus2 nunu nurh rhs2 rhnu rhrh]
            f_seq=[f1(:,1).*f2 f1(:,2).*f2 f1(:,3).*f2];
            s_1=s_1+sum(cdd_l2.*f_seq./seq_ep,1);

            % Arthur tends to comment this out and double s_1 instead. Are s_1
            % and s_2 always the same?
            f1=f(ep_ind3,:);
            f_seq2=[f1(:,1).*flipud(f1) f1(:,2).*flipud(f1) f1(:,3).*flipud(f1)];
            s_2=s_2+sum(cad_l2.*f_seq2./seq_ep2,1);
          end
        end
      else 
       % Useful for debugging; calculation time depends on grid-size
       % Display progress?
       dispprog=0;
       for my=(-Ny+1:Ny-1)
         if dispprog
           disp(sprintf('%i%% complete',...
                floor(((my+Ny-1)/(2*Ny)*100))))
         end
         for mx=(-Nx+1:Nx-1)                                                    
           % Create two logical vectors for indexing the offset in diagonals
           % of the FFT of the autocovariance
           a=y<Ny-my&x<Nx-mx;                                                
           b=y>=-my&x>=-mx;                                                  
           ep_ind1=a&b;                                                      
           ep_ind1=ep_ind1(:);
           a=y>=my&x>=mx;                                                    
           b=y<Ny+my&x<Nx+mx;                                                
           ep_ind2=a&b;                                                      
           ep_ind2=ep_ind2(:);

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
           fftcov_diag=blurosy(th,paramsSbar,xver,'efd',[my mx]);                        
           fftcov_diag=v2s(fftcov_diag,params)*nf;                                           
           cdd=fftcov_diag(max([1,my+1]):min([Ny+my,Ny]),...                         
                           max([1,mx+1]):min([Nx+mx,Nx]));                           
           cdd_l2=abs(cdd).^2;                                                  
           cdd_l2=cdd_l2(:);                                                       
                                                                                 
           % Use the same framework to select the antidiagonals from BLUROSY
           fftcov_anti=blurosy(th,paramsSbar,xver,'efd',[my2 mx2]);
           fftcov_anti=v2s(fftcov_anti,params)*nf;
           cad=fftcov_anti(max([1,my2-Ny+2]):min([my2+1,Ny]),...
                           max([1,mx2-Nx+2]):min([mx2+1,Nx]));
           cad_l2=abs(cad).^2;
           cad_l2=cad_l2(:);

           if normalize
             seq_ep =Sbar_n(ep_ind1).*Sbar_n(ep_ind2);
             ep_n3=ep_n(ep_ind3);
             seq_ep2=ep_n3.*flip(ep_n3,1);
           end

           % The partials indexed for the current `diagonal'
           f1=f(ep_ind1,:);f2=f(ep_ind2,:);
           % Product of all partials in row-order:
           % [s2s2 s2nu s2rh nus2 nunu nurh rhs2 rhnu rhrh]
           f_seq=[f1(:,1).*f2 f1(:,2).*f2 f1(:,3).*f2];
           s_1=s_1+sum(cdd_l2.*f_seq./seq_ep,1);

           % Arthur tends to comment this out and double s_1 instead. Are s_1
           % and s_2 always the same?
           f1=f(ep_ind3,:);
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
     % Following several tests, we see that there is a very small difference 
     % between s_1 and s_2 following summation (O(-10) - O(-16)), with 
     % dependence on th.
     % disp(sprintf('s_1 - s_2: %g %g %g %g %g %g %g %g %g',s_1-s_2))
     s=reshape(s_1+s_2,np,np);
     if isfield(params,'taper')&numel(params.taper)==prod(NyNx)
        normfact=(prod(params.NyNx)./sum(Tx(:).^2)).^2.*prod(params.NyNx).^2;
     else
        normfact=prod(params.NyNx).^2;
     end
     covg=s./normfact;
     % wkspc=whos;sum([wkspc.bytes])
     grads=NaN;
    end
  covg=ifinvslc(covg,ifinv);

  % Optional output
  varns={covg,grads};
  varargout=varns(1:nargout);
elseif strcmp(th,'demo1')
    % Makes some plots of the covariance, periodogram, and their partial 
    % derivatives  wrt the Matern parameters
    th=[1 0.5 1];params=[];params.NyNx=[32 32];params.dydx=[1 1];
    if ismember(th(2),[1/3,1/2,1,3/2,5/2,Inf])
       [dSbards2,~,~,d_acv_s2]=blurosy(th,params,[],[],[],1);
       [dSbardnu,~,~,d_acv_nu]=blurosy(th,params,[],[],[],2);
       [dSbardrh,~,~,d_acv_rh]=blurosy(th,params,[],[],[],3);
       clf;subplot(331);imagefnan(v2s(p));title('|Hk^2|')
       subplot(334);imagefnan(v2s(Sbar));title('Sk')
       subplot(337);imagefnan(v2s(acv));title('Cy')
       subplot(332);imagefnan(v2s(dSbards2));title('dSkds2');
       subplot(333);imagefnan(v2s(d_acv_s2));title('dCyds2');
       subplot(335);imagefnan(v2s(dSbardnu));title('dSkdnu');
       subplot(336);imagefnan(v2s(d_acv_nu));title('dCydnu');
       subplot(338);imagefnan(v2s(dSbardrh));title('dSkdrh');
       subplot(339);imagefnan(v2s(d_acv_rh));title('dCydrh');
     else
       error('we are not ready for the general Matern case yet, choose a special value of nu')
     end
elseif strcmp(th,'demo2')
    % Figures developed for direct comparison with Arthur's covgammiosl_sample
    % Need to acquire an MLEOSL suite
    % th1=[1 1/2 5]; 
    %th1=[10 1/2 10]; 
    th1=[5 1/2 7]; 
    params=[];params.NyNx=[94 97];params.dydx=[1 1];
    params.blurs=Inf;params.taper=1;
    ifinv=[1 0 1];
    try 
      [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(date,100);
    catch
      try
        mleosl('demo1',0,th1,params,[],[],th1,[1 0 1])
      catch
        mleosl('demo1',200,th1,params,[],[],th1,[1 0 1])
      end
      [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(date,100);
    end
    covemp=cov(thhats);
    mobs=nanmean(thhats);

    % For the covgammiosl sampling method, we want to study the effect of increasing
    % number of samples and will recalculate the variance using F1inv depending
    % on the number of samples
    [~,~,~,k]=simulosl(th1,params,1);
    F1=fishiosl(k,th1,0,params);
    F1=ifinvslc(F1,ifinv);
    F1inv=inv(F1);
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
        covg1=ifinvslc(covg1,[1 0 1]);
        cov1=F1inv*covg1*F1inv;
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

    % Plot it up
    figure(1);
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
            'LineWidth',2,'DisplayName','empirical');
    nspl(1)=plot(sigvec,normpdf(sigvec,th1(:,1),sqrt(cov1s2s2(end))),...
            'LineWidth',2,'DisplayName','sample');
    ndpl(1)=plot(sigvec,normpdf(sigvec,th1(:,1),sqrt(cov2(1,1))),':',...
            'LineWidth',2,'DisplayName','diagonal');
    nfpl(1)=plot(sigvec,normpdf(sigvec,th1(:,1),sqrt(cov3(1,1))),'--',...
            'LineWidth',2,'DisplayName','dftmtx');
    xlabel('\sigma^2')
    ylabel('density')
    legend()

    % Histogram of rh estimates
    subplot(222)
    hpl(2)=histogram(thhats(:,3),numbins,'Normalization','pdf');
    vpl(2)=xline(th1(:,3),'Label','\theta_0','LabelOrientation','horizontal');
    vecmod=[0.95 1.05];
    veclen=numsims*10;
    hold on
    nepl(2)=plot(rhovec,normpdf(rhovec,mobs(:,3),sqrt(covemp(3,3))),...
                'LineWidth',2);
    nspl(2)=plot(rhovec,normpdf(rhovec,th1(:,3),sqrt(cov1rhrh(end))),...
                'LineWidth',2);
    ndpl(2)=plot(rhovec,normpdf(rhovec,th1(:,3),sqrt(cov2(end,end))),':',...
                'LineWidth',2);
    nfpl(2)=plot(rhovec,normpdf(rhovec,th1(:,3),sqrt(cov3(end,end))),'--',...
                'LineWidth',2);
    xlabel('\rho')


    % Calculations for error ellipses
    clrs=colororder;
    t=linspace(0,2*pi);
    cl=0.95;
    s=chi2inv(cl,2);
    [V,D]=eig(ifinvslc(covemp,ifinv));
    aemp=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    [V,D]=eig(cov1);
    asample=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    [V,D]=eig(cov2);
    a2=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    [V,D]=eig(cov3);
    a3=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
    % Cross plot with error ellipses
    subplot(223)
    spl(1)=plot(thhats(:,1),thhats(:,3),'o');
    hold on
    er(1)=plot(aemp(1,:)+mobs(1),aemp(2,:)+mobs(3),'Color',clrs(2,:));
    er(2)=plot(asample(1,:)+mobs(1),asample(2,:)+mobs(3),'Color',clrs(3,:));
    er(3)=plot(a2(1,:)+mobs(1),a2(2,:)+mobs(3),'Color',clrs(4,:),':');
    er(4)=plot(a3(1,:)+mobs(1),a3(2,:)+mobs(3),'Color',clrs(5,:),'--');
    vpl(3)=xline(th1(:,1),'Label','\sigma^2_0','LabelOrientation','horizontal');
    hpl(3)=yline(th1(:,3),'Label','\rho_0','LabelOrientation','horizontal');
    xlabel('$\hat{\sigma^2}$','Interpreter','latex')
    ylabel('$\hat{\rho}$','Interpreter','latex')
    grid on
    box on

    % Covariance comparison
    subplot(224)
    spl(2)=plot(1:numsims,cov1s2s2,'-','Color',clrs(3,:),'DisplayName','{s2,s2}');
    hold on
    spl(3)=plot(1:numsims,cov1s2rh,'--','Color',clrs(3,:),'DisplayName','{s2,rh}');
    spl(4)=plot(1:numsims,cov1rhrh,':','Color',clrs(3,:),'DisplayName','{rh,rh}');
    hJ(1)=yline(covemp(1,1),'-','Color',clrs(2,:));
    hJ(2)=yline(covemp(1,end),'--','Color',clrs(2,:));
    hJ(3)=yline(covemp(end,end),':','Color',clrs(2,:));
    hJ(4)=yline(cov2(1,1),'-','Color',clrs(4,:));
    hJ(5)=yline(cov2(1,end),'--','Color',clrs(4,:));
    hJ(6)=yline(cov2(end,end),':','Color',clrs(4,:));
    hJ(7)=yline(cov3(1,1),'-','Color',clrs(5,:));
    hJ(8)=yline(cov3(1,end),'--','Color',clrs(5,:));
    hJ(9)=yline(cov3(end,end),':','Color',clrs(5,:));
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
      sprintf('covgammiosl_demo2_parametercovariancemethods_%ix%i_s2%irh%i.eps',...
        params.NyNx.*params.dydx,th1(1),th1(3)),'epsc')

    % Normalized covariance heatmaps
    figure(2)
    kelicol
    cmrange=[0.5 1];
    labs={'s2','nu','rho'};
    labs=ifinvslc(labs,ifinv);
    np=sum(ifinv);
    ah(1)=subplot(221);
    im(1)=imagesc(ifinvslc(covemp,ifinv)./...
          [diag(sqrt(ifinvslc(covemp,ifinv)))*...
           diag(sqrt(ifinvslc(covemp,ifinv)))'],...
           cmrange);
    title('empirical')
    axis image

    ah(2)=subplot(222);
    im(2)=imagesc(cov1./[diag(sqrt(cov1))*diag(sqrt(cov1))'],cmrange);
    title('sample')
    axis image

    ah(3)=subplot(223);
    im(3)=imagesc(cov2./[diag(sqrt(cov2))*diag(sqrt(cov2))'],cmrange);
    title('diagonal')
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

elseif strcmp(th,'demo3')
  % Sanity check developed for evaluating the calculation of the periodogram's
  % covariance in method 3; 
  %%% OLW TODO: NEEDS TO BE UPDATED

  % Push back one argument
  th=params;
  grd(1).dydx=[1 1];grd(2).dydx=[1 1];
  grd(1).NyNx=[127 128];grd(2).NyNx=[128 128];

  for ind=1:3
    if ind<3
      % The distance grid, dxxp
      dydx=grd(ind).dydx;NyNx=grd(ind).NyNx;
      if NyNx(1)==NyNx(2)
        ydim=[-floor(NyNx(1)/2):ceil(NyNx(1)/2)-1]';
        xdim=ydim;
        dxxp{ind}=xxpdist([xdim ydim]); 
      else
        ydim=[-floor(NyNx(1)/2):ceil(NyNx(1)/2)-1]';
        xdim=[-floor(NyNx(2)/2):ceil(NyNx(2)/2)-1]';
        ydim=[0:NyNx(1)-1]';
        xdim=[0:NyNx(2)-1]';
        keyboard

        [MDX,NDXP]=meshgrid(xdim(:),ydim(:));
        %dxxpf=xxpdist([MDX(:) NDXP(:)],[MDX(1) NDXP(1)]);
        dxxpf=xxpdist(MDX(:), NDXP(:));
        dxxp{ind}=v2s(dxxpf,grd(ind));
      end
      % The covariance matrix 
      Cmn{ind}=maternosy(dxxp{ind},th);
      % The DFT matrix (mxm)
      Um{ind}=dftmtx(size(Cmn{ind},2)); 
      % The DFT matrix (nxn)
      Un{ind}=dftmtx(size(Cmn{ind},1)); 
      % The covariance of the periodogram via Isserlis' rule
      Hk2covf=abs(Um{ind}'*Cmn{ind}*Un{ind}).^2+...
        abs(Um{ind}'*Cmn{ind}*conj(Un{ind})).^2;
      Hk2cov{ind}=Hk2covf/prod(NyNx.*dydx);
    else
      NyNx3=[min([grd(1).NyNx(1) grd(2).NyNx(1)]),...
             min([grd(1).NyNx(2) grd(2).NyNx(2)])];
      dxxp{ind}=dxxp{1}(1:NyNx3(1),1:NyNx3(2))-dxxp{2}(1:NyNx3(1),1:NyNx3(2));
      Cmn{ind}=Cmn{1}(1:NyNx3(1),1:NyNx3(2))-Cmn{2}(1:NyNx3(1),1:NyNx3(2));
      Hk2cov{ind}=Hk2cov{1}(1:NyNx3(1),1:NyNx3(2))-Hkcov{2}(1:NyNx3(1),1:NyNx3(2));
    end

    figure(ind);
    clf
    ah(1)=subplot(221);imagesc(dxxp{ind});axis image ij;
    if ind==1
      cl{1}=clim;
    else
      set(gca,'CLim',cl{1});
    end
    title('2-D lag-distance grid')
    xlabel('x1 index')
    ylabel('x2 index')
    ah(2)=subplot(222);imagesc(Cmn{ind});axis image ij;
    if ind==1
      cl{2}=clim;
    else
      set(gca,'CLim',cl{2});
    end
    title('Spatial covariance, Cx')
    xlabel('x1 index')
    ylabel('x2 index')
    ah(3)=subplot(223);imagesc(Hkcov{ind});axis image ij;
    if ind==1
      cl{3}=clim;
    else
      set(gca,'CLim',cl{3});
    end
    title('cov(I_n(\omega_1),I_n(\omega_2))')
    xlabel('\omega_1 index')
    ylabel('\omega_2 index')
    if ind<3
        ah(4)=subplot(224);imagesc(log10(Hkcov{ind}));axis image ij;
        title('log10(cov(I_n(\omega_1),I_n(\omega_2)))')
    else
        ah(4)=subplot(224);imagesc(log10(abs(Hkcov{ind})));axis image ij;
        title('log10(abs(cov(I_n(\omega_1),I_n(\omega_2))))')
    end
    if ind==1
      cl{4}=clim;
    else
      set(gca,'CLim',cl{4});
    end
    xlabel('\omega_1 index')
    ylabel('\omega_2 index')
    if ind<3
      ti{ind}=sgtitle(sprintf('%ix%i grid | th=[%0.2f %0.2f %0.2f]',NyNx,th));
    else
      ti{ind}=sgtitle({sprintf(...
         'Difference between %ix%i and %ix%i grid calculations',...
           grd(1).NyNx,grd(2).NyNx),
         sprintf('th=[%0.2f %0.2f %0.2f]',th)});
    end
    if ind<3
      fnam=sprintf('covperiodogram_%iby%i_%0.0f_%0.0f_%0.0f',NyNx,th);
    else
      fnam=sprintf('covperiodogram_griddiff_%0.0f_%0.0f_%0.0f',th);
    end
    print('-depsc',fnam)
  end
elseif strcmp(th,'demo4')
    % Simulate fields, estimate Matern parameters, and calculate their
    % covariance, all as a function of the NUMBER OF SIMULATIONS; this is
    % mostly just a validation of the COVGAMMIOSL methods

    % Assign the true Matern parameters that we will simulate from on our way to
    % forming the empirical covariance. We will invert for all three parameters.
    th0=[5 0.5 5];
    ifinv=[1 1 1];
    % Create a grid that is large enough for positive definite embedding in SGP,
    % but small enough that we can form the entire DFTMTX (i.e., less than about
    % 1e4 grid points). We will assume complete and evenly sampled field.
    params=[];params.NyNx=[58 53];params.dydx=[1 1];                            
    params.blurs=Inf;params.taper=1;                                            
    % NOTE: Try varying dydx in the future to diagnose whether this has been
    % incorporated through the COVGAMMIOSL methods correctly

    % Begin the simulation study first using MLEOSL's DEMO1 to create a suite of
    % realizations from th0 for the grid defined by params


    % Calculate the empirical mean and covariance
    mobs  =mean(thhats); 
    covobs=cov(thhats);

    % Now, decide whether to calculate the analytical covariance for th0, the
    % empirical mean, or ANY given thhat, which is more true to application
    thcov=th0;%mobs;%thhats(randi(size(thhats,1)),:);
     
    % Set up the 'SAMPLING' method (1), which outputs the gradient of the score
    % GRADS for every data vector sampled for THCOV, allowing us to see the
    % evolution of the covariance approximation with increasing numbers of
    % samples. Note that these samples are obviously independent from the 
    % realizations of the MLEOSL DEMO1 suite.
    [covg1,grads]=covgammiosl(th1,params,1,ifinv);
    [~,~,~,k]=simulosl(th1,params,1);
    F1=fishiosl(k,th1,0,params);
    F1=ifinvslc(F1,ifinv);
    F1inv=inv(F1);
    numsims=size(grads,1);
    for ind=2:numsims
        covg1=real(cov(grads(1:ind,:)))./prod(params.NyNx).^2;
        covg1=ifinvslc(covg1,[ifinv]);
        cov1=F1inv*covg1*F1inv;
        cov1(:,:,ind)=trilos(cov1);
    end
    %%% OLW TODO: NEEDS TO BE UPDATED

    covg2  =covgammiosl(th1,params,2,ifinv);
    covg3  =covgammiosl(th1,params,3,ifinv);
    cov2=varianceofestimates(th1,params,covg2,ifinv);
    cov2=trilos(cov2);
    cov3=varianceofestimates(th1,params,covg3,ifinv);
    cov3=trilos(cov3);
    labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
          '\sigma^2,\rho', '\nu,\rho'};

    figure(1)
    [ah,ha,H]=krijetem(subnum(2,3));

    for ind=2:25:numsims
        covg1=real(cov(grads(1:ind,:)))./prod(params.NyNx).^2;
        covg1=ifinvslc(covg1,[ifinv]);
        cov1=F1inv*covg1*F1inv;
        cov1=trilos(cov1);
        for jnd=1:numel(cov1)
          axes(ah(jnd))
          plot(ind,cov1(jnd),'ko')
          hold on
        end
    end
    for jnd=1:numel(cov1)
      axes(ah(jnd))
      if jnd==1|jnd==4
          ylabel('Covariance')
      end
      if jnd>3
          xlabel('Number of simulations')
      end
      yline(cov2(jnd),'DisplayName','diagonals')
      yline(cov3(jnd),'DisplayName','dftmtx')
      title(sprintf('%s',labs{jnd}))
    end
   ti=sgtitle(sprintf('\\sigma^2=%0.2f, \\nu=%0.2f, \\rho=%0.2f; %im x %im',...
           th1,params.NyNx.*params.dydx))
   movev(ah,-0.03)

   saveas(gcf,...
     sprintf('covgammiosl_demo4_1_paramcovnumsims_%ix%i_s2%inu%irh%i.eps',...
       params.NyNx.*params.dydx,th1),'epsc')

elseif strcmp(th,'demo5')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the GRID SIZE, which is presented in terms
   % of pi rho by pi rho. Plot the bias from the empirical study for each
   % parameter.
   params=[];params.dydx=[1 1];params.blurs=Inf;params.taper=1;ifinv=[1 1 1];
   numempsims=200;

   % Labels
   labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
          '\sigma^2,\rho', '\nu,\rho'};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];

   figure(2)
   clf;
   [ah,ha,H]=krijetem(subnum(3,3));
   % For now, small values so that we have a chance at forming the data vector
   % and finding an estimate for small grids 
   rh=2;
   th1=[1 0.5 rh];
   % Sample in logspace
   loggs=logspace(0,1,10);
   thhats=zeros(numempsims,3,size(loggs,2));
   try
     % If we have already calculate the estimates, load them. Note that we may
     % need to append a date to the file name
     datum='13-Mar-2025';%date;
     thhats=load(sprintf('covgammiosl_demo5_empthhats_%s.mat',datum),'thhats');
     thhats=thhats.thhats;
     % Consider trimming; some of results are truly bad on small grids
     % e.g., [gthhats,trimi]=trimit(thhats(:,:,1),80);
   catch
     for ind=1:size(loggs,2)
       % Create a larger grid with length linearly increasing in pi rho
       gs=rh*pi*loggs(ind); 
       params.NyNx=[floor(gs) floor(gs)];
       for mnd=1:numempsims
         Hx=simulosl(th1,params,0);
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
       save(sprintf('covgammiosl_demo5_empthhats_%s.mat',date),'thhats')
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
     cov1c=varianceofestimates(th1,params,covg1,ifinv);
     % Calculate the 'dftmtx' covariance
     covg2  =covgammiosl(th1,params,2,ifinv);
     cov2c=varianceofestimates(th1,params,covg2,ifinv);
     % Calculate the 'diagonals' covariance
     covg3  =covgammiosl(th1,params,3,ifinv);
     cov3c=varianceofestimates(th1,params,covg3,ifinv);
     try
       % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
       covempt(ind,:)=trilos(covemp);
       cov1t(ind,:)  =trilos(cov1c);
       cov2t(ind,:)  =trilos(cov2c);
       cov3t(ind,:)  =trilos(cov3c);
     catch
       keyboard
     end
   end
   cmnylim=[minmax([abs(covempt) abs(cov1t) abs(cov2t) abs(cov3t)])];
   cmnylimrng=round(log10(cmnylim),0);
   for jnd=1:6
     axes(ah(jnd))
     % Plot reference diagonals as light gray lines on bottom for each power of
     % 10 between cmnylim
     hold on
     keyboard
     rnd=0;
     for ind=cmnylimrng(1):cmnylimrng(2)
       rnd=rnd+1;
       refd(rnd)=loglog(loggs,10^(ind)./floor(rh.*pi.*loggs) ,'Color',[0.8 0.8 0.8]);
     end
     % Plot the absolute covariance so that negative correlations do not 
     % mislead
     loglog(loggs,abs(covempt(:,jnd)),'ko','DisplayName','empirical')
     loglog(loggs,abs(cov1t(:,jnd)),'^','MarkerEdgeColor',clrs{1},...
            'MarkerFaceColor',clrs{1},'DisplayName','sample')
     loglog(loggs,abs(cov2t(:,jnd)),'s','MarkerEdgeColor',clrs{2},...
            'MarkerFaceColor',clrs{2},'DisplayName','dftmtx')
     loglog(loggs,abs(cov3t(:,jnd)),'*','MarkerFaceColor',clrs{3},...
            'MarkerFaceColor',clrs{3},'DisplayName','diagonal')
   end
   for jnd=1:6
     axes(ah(jnd));
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
     loglog(loggs,10.^(predemp(jnd,:)),'Color','k');
     loglog(loggs,10.^(predcov1(jnd,:)),'Color',clrs{1});
     loglog(loggs,10.^(predcov2(jnd,:)),'Color',clrs{2});
     loglog(loggs,10.^(predcov3(jnd,:)),'Color',clrs{3});

     % Annotate the slopes
     % Not enough space for the following:
     % leg1(ind)=legend(sprintf('emp, %s%0.1f}','(\pi\rho)^{',mbfemp(1,1)),...
     %                  sprintf('sam, %s%0.1f}','(\pi\rho)^{',mbfcov1(1,1)),...
     %                  sprintf('dft, %s%0.1f}','(\pi\rho)^{',mbfcov2(1,1)),...
     %                  sprintf('dia, %s%0.1f}','(\pi\rho)^{',mbfcov3(1,1)),...
     %                  'interpreter','tex','BackgroundAlpha',0,'box','off');
     % instead, just report the slope and describe in a figure caption:
     keyboard
     [leg1(jnd),legic]=legend('','','',...
                      sprintf('empirical, %0.1f',abs(mbfemp(jnd,1))),...
                      sprintf('sample, %0.1f',abs(mbfcov1(jnd,1))),...
                      sprintf('dftmtx, %0.1f',abs(mbfcov2(jnd,1))),...
                      sprintf('diagonal, %0.1f',abs(mbfcov3(jnd,1))),...
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
     title(sprintf('%s',labs{jnd}))
     longticks
     ylim([cmnylim(1)*0.5 cmnylim(2)*2.0])
   end      
    
   % Plot the bias of the estimates from the empirical study
   axes(ah(7));
   semilogx(loggs,squeeze(mean(thhats(:,1,:))-th1(1)),'ko')
   % Add the 0-line where the estimate equals the truth
   yline(0)
   xlabel('Length ($\pi\rho$)','Interpreter','latex')
   ylabel('$\langle\hat{\theta}\rangle-\theta_0$','interpreter','latex');
   ti(7)=title('\sigma^2','FontWeight','bold','Interpreter','tex');
   longticks
   %
   axes(ah(8));
   semilogx(loggs,squeeze(mean(thhats(:,2,:))-th1(2)),'ko')
   yline(0)
   xlabel('Length ($\pi\rho$)','Interpreter','latex')
   ti(8)=title('\nu','FontWeight','bold','Interpreter','tex');
   ylim([-0.1,max(mean(thhats(:,2,:))-th1(2))+0.1])
   longticks
   %
   axes(ah(9));
   semilogx(loggs,squeeze(mean(thhats(:,3,:))-th1(3)),'ko')
   yline(0)
   xlabel('Length ($\pi\rho$)','Interpreter','latex')
   ti(9)=title('\rho','FontWeight','bold','Interpreter','tex');
   longticks
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
     sprintf('covgammiosl_demo5_paramcovgridsize_s2%inu%irh%i_%s.eps',...
       th1.*[1 10 1],date),'epsc')
elseif strcmp(th,'demo6')
% Simulate fields, estimate Matern parameters, and calculate their
% covariance, all as a function of the percent of RANDOM DELETIONS
th1=[1 0.5 2];
params.NyNx=[36 36];params.dydx=[1 1];params.blurs=Inf;
params.blurs=Inf;params.taper=1;                                            
oneprcnt=ceil(0.01.*prod(params.NyNx))
ifinv=[1 1 1];                                                              

labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
      '\sigma^2,\rho', '\nu,\rho'};
% color friendly palette
clrs=["#1B9E77","#D95F02","#7570B3"];

clf;
taper=ones(params.NyNx);
[ah,ha,H]=krijetem(subnum(2,3));
%for ind=1:20
for ind=1:10
 thhats=[];
 params.taper=taper;
 prcmiss=1-sum(taper,"all")./prod(params.NyNx);
 for mnd=1:100
    Hx=simulosl(th1,params,0);
    thhat=NaN;
    while isnan(thhat)
      [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
    end
    thhats(mnd,:)=thhat.*scl;
 end
 keyboard
 covemp=nancov(thhats);
 covg1  =covgammiosl(th1,params,1,ifinv);
 covg2  =covgammiosl(th1,params,2,ifinv);
 covg3  =covgammiosl(th1,params,3,ifinv);
 cov1c=varianceofestimates(th1,params,covg1,ifinv);
 cov2c=varianceofestimates(th1,params,covg2,ifinv);
 cov3c=varianceofestimates(th1,params,covg3,ifinv);
 try
   cov1t  =trilos(cov1c);
   cov2t  =trilos(cov2c);
   cov3t  =trilos(cov3c);
   covempt=trilos(covemp);
 catch
   keyboard
 end
 for jnd=1:6
    axes(ah(jnd))
    plot(prcmiss*100,covempt(jnd),'ko','DisplayName','empirical')
    hold on
    plot(prcmiss*100,cov1t(jnd),'^','MarkerEdgeColor',clrs(1),...
         'MarkerFaceColor',clrs(1),'DisplayName','sample')
    plot(prcmiss*100,cov2t(jnd),'s','MarkerEdgeColor',clrs(2),...
         'MarkerFaceColor',clrs(2),'DisplayName','dftmtx')
    plot(prcmiss*100,cov3t(jnd),'*','MarkerEdgeColor',clrs(3),...
         'MarkerFaceColor',clrs(3),'DisplayName','diagonal')
 end
 %taper(randi(prod(params.NyNx),1,50))=deal(0);
 taper(randi(prod(params.NyNx),1,oneprcnt))=deal(0);
end
for jnd=1:6
 axes(ah(jnd))
 if jnd==1|jnd==4
    ylabel('$\mathrm{cov}\{\theta_i,\theta_j\}$','Interpreter','latex')
 end
 if jnd>3
    xlabel('Percent random deletion')
 end
 title(sprintf('%s',labs{jnd}))
end      
keyboard
ti=sgtitle(...
sprintf('$%s%0.2f%s%0.2f%s%0.2f$%s%i%s%i%s',...
'\sigma^2=',th1(1),' \mathrm{[unit]}^2, \nu=',th1(2),', \rho=',th1(3),...
' m; ',params.NyNx(1).*params.dydx(1),' m x ',...
params.NyNx(2).*params.dydx(2),'m'),... 
       'Interpreter','latex');
movev(ah,-0.03)                                                              
                                                                            
saveas(gcf,...                                                               
 sprintf('covgammiosl_demo6_samp_paramcovprcmiss_%i_%i_%i_%i_%i_%s.eps',...       
   th1.*[1 10 1],params.NyNx.*params.dydx,date),'epsc')  
elseif strcmp(th,'demo7')
   keyboard
   % mleosl('demo2','18-Feb-2025')
   % [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload('18-Feb-2025',100);

   p=[];p.NyNx=[91 91];p.dydx=[1 1];p.blurs=Inf;p.taper=1;
   th0=[3 1.5 5];
   % % mleosl('demo1',300,[3 1.5 5],p)
   % %%% OLW -- just load it straight from mleosl_thhat_05-Mar-2025
   thhats=load('mleosl_thhat_12-Mar-2025');
keyboard
   % color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];

   % p=[];p.NyNx=[64 64];p.dydx=[1 1];p.blurs=Inf;p.taper=1;
   % th0=[2 0.8 3];
   % mleosl('demo1',300,th0,p,[],[],th0,[1 1 1])
   %%% OLW -- just load it straight from mleosl_thhat_05-Mar-2025
   % thhats=load('mleosl_thhat_05-Mar-2025');
   % thhats=thhats(264:end,:);

   covg1=covgammiosl(th0,p,1,[1 1 1])
   cov1=varianceofestimates(th0,p,covg1,[1 1 1])
   cov1s2nu=ifinvslc(cov1,[1 1 0]);
   cov1s2rh=ifinvslc(cov1,[1 0 1]);
   cov1nurh=ifinvslc(cov1,[0 1 1]);

   covg2=covgammiosl(th0,p,2,[1 1 1])
   cov2=varianceofestimates(th0,p,covg2,[1 1 1])
   cov2s2nu=ifinvslc(cov2,[1 1 0]);
   cov2s2rh=ifinvslc(cov2,[1 0 1]);
   cov2nurh=ifinvslc(cov2,[0 1 1]);

   covg3=covgammiosl(th0,p,3,[1 1 1])
   cov3=varianceofestimates(th0,p,covg3,[1 1 1])
   cov3s2nu=ifinvslc(cov3,[1 1 0]);
   cov3s2rh=ifinvslc(cov3,[1 0 1]);
   cov3nurh=ifinvslc(cov3,[0 1 1]);

   cove=cov(thhats);
   coves2nu=ifinvslc(cove,[1 1 0]);
   coves2rh=ifinvslc(cove,[1 0 1]);
   covenurh=ifinvslc(cove,[0 1 1]);

   mobs1=nanmean(thhats(:,[1 2]))
   mobs2=nanmean(thhats(:,[1 3]))
   mobs3=nanmean(thhats(:,[2 3]))

   t=linspace(0,2*pi);
   cl=0.95;
   s=chi2inv(cl,2);

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
   % need first axis
   subplot(131)
   plot(thhats(:,1),thhats(:,2),'ko')
   hold on
   keyboard
   e1(1)=plot(a11(1,:)+mobs1(1),a11(2,:)+mobs1(2),'Color',clrs(1),'LineWidth',1);
   e2(1)=plot(a21(1,:)+mobs1(1),a21(2,:)+mobs1(2),'Color',clrs(2),'LineWidth',1);
   e3(1)=plot(a31(1,:)+mobs1(1),a31(2,:)+mobs1(2),'Color',clrs(3),'LineWidth',1);
   ee(1)=plot(ae1(1,:)+mobs1(1),ae1(2,:)+mobs1(2),'Color','k','LineWidth',1);
   longticks
   xlabel('variance, $\sigma^2$','Interpreter','latex')
   ylabel('smoothness, $\nu$','Interpreter','latex')
   % need second axis
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
   % need third axis
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
           th0,p.NyNx.*p.dydx),'interpreter','latex')

   keyboard
   saveas(gcf,...
     sprintf('covgammiosl_demo7_crossplots_%i_%i_%i_%i_%i_%s.eps',...
       th0.*[1 10 1],p.NyNx.*p.dydx,date),'epsc')
   
   % Normalized covariance heatmaps
   figure(20)
   kelicol
   cmrange=[-1 1];
   %cmrange=[0.5 1];
   labs={'s2','nu','rho'};
   labs=ifinvslc(labs,ifinv);
   np=sum(ifinv);
   ah(1)=subplot(221);
   im(1)=imagesc(ifinvslc(cove,ifinv)./...
         [diag(sqrt(ifinvslc(cove,ifinv)))*...
          diag(sqrt(ifinvslc(cove,ifinv)))'],...
          cmrange);
   title('empirical')
   axis image

   ah(2)=subplot(222);
   im(2)=imagesc(cov1./[diag(sqrt(cov1))*diag(sqrt(cov1))'],cmrange);
   title('sample')
   axis image

   ah(3)=subplot(223);
   im(3)=imagesc(cov2./[diag(sqrt(cov2))*diag(sqrt(cov2))'],cmrange);
   title('diagonal')
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

   moveh(ah(2),-0.04)
   movev(ah(1),-0.02)
   ah(2).Position(4)=ah(1).Position(4);ah(2).Position(3)=ah(1).Position(3);
   ah(2).Position(2)=ah(1).Position(2);
   ti=sgtitle(sprintf('\\sigma^2=%0.2f, \\nu=%0.2f, \\rho=%0.2f; %im x %im',...
           th0,p.NyNx.*p.dydx))

   saveas(gcf,...
     sprintf('covgammiosl_demo7_corrplots_%ix%i_s2%inu%irh%i.eps',...
       p.NyNx.*p.dydx,th0),'epsc')
elseif strcmp(th,'demo8')
   % Simulate fields, estimate Matern parameters, and calculate their
   % covariance, all as a function of the GRID SPACING to study INFILL
   % ASYMPTOTICS
   params=[];params.blurs=Inf;params.taper=1;ifinv=[1 1 1];
   numempsims=3;

   % Labels
   labs={'\sigma^2,\sigma^2', '\nu,\nu', '\rho,\rho', '\sigma^2,\nu',...
          '\sigma^2,\rho', '\nu,\rho'};
   % Color friendly palette
   clrs=["#1B9E77","#D95F02","#7570B3"];

   figure(2)
   clf;
   [ah,ha,H]=krijetem(subnum(3,3));
   % Set the Matern parameters 
   rh=2;
   th1=[1 0.5 rh];
   % Case 1:
   % Do we want to have the same number of samples in each experiment and to be
   % able to calculate integer NyNx from increasing integer dydx? For this to be
   % possible, N^2 must be factorable by all dydx values.
   % A candidate constant grid size is 40, just less than 7 correlation lengths
   % for rho 2, and dydxs=[1 2 4 8 10].
   % Case 2:
   % If we are increasing the sample size of the grid with each increase in
   % dydx (i.e., constant NyNx, not constant N^2), then we might as well sample
   % dydxs in logspace to have a visually consistent plot with demo5, even
   % though we are not expecting to see a linear trend in loglog here (nor much
   % of any trend). However, this is increasing the number of samples similar to
   % demo5, so maybe we do want the above case...
   cs=1;
   if cs==1
     % Constant grid size (dydx*NyNx=constant) for all experiments 
     gslen=40;
     % Assuming a square grid and providing a single length for each experiment
     dydxs=[1 2 4 8 10];
     NyNxs=gslen./dydxs;
   elseif cs==2 
     % Constant number of samples (NyNx=constant) for all experiments 
     % Calculate increasing sample size in logspace
     loggs=logspace(0,1,5);
     % Assuming a square grid and providing a single length for each experiment
     dydxs=ceil(loggs);
     NyNxs=repelem(ceil(pi*rh*7),size(dydxs,1)); 
     gslen=dydxs.*ceil(pi*rh*7);
   end
   % Calculate the ensemble
   thhats=zeros(numempsims,3,size(dydxs,2));
   try
     % If we have already calculate the estimates, load them. Note that we may
     % need to append a date to the file name
     thhats=load(sprintf('covgammiosl_demo8_empthhats_%s.mat',date),'thhats');
     thhats=thhats.thhats;
   catch
     for ind=1:size(dydxs,2)
       % Set spacing and number of samples
       params.NyNx=[NyNxs(ind) NyNxs(ind)];
       params.dydx=[dydxs(ind) dydxs(ind)];
       for mnd=1:numempsims
         Hx=simulosl(th1,params,0);
         thhat=NaN;
         while isnan(thhat)
           [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],th1,ifinv,0);
         end
         thhats(mnd,:,ind)=thhat.*scl;
       end
     end
     % Consider saving your calculations by setting 'svthhats' in the debugger
     keyboard
     svthhats=0;
     if svthhats
       save(sprintf('covgammiosl_demo8_empthhats_%s.mat',date),'thhats')
     end
   end
   % Form the covariance matrices
   for ind=1:size(dydxs,2)
     % Set spacing and number of samples
     params.NyNx=[NyNxs(ind) NyNxs(ind)];
     params.dydx=[dydxs(ind) dydxs(ind)];
     % Calculate the empirical covariance from the ensemble
     covemp=nancov(thhats(:,:,ind));
     % Calculate the 'sample' covariance
     covg1  =covgammiosl(th1,params,1,ifinv);
     cov1c=varianceofestimates(th1,params,covg1,ifinv);
     % Calculate the 'dftmtx' covariance
     covg2  =covgammiosl(th1,params,2,ifinv);
     cov2c=varianceofestimates(th1,params,covg2,ifinv);
     % Calculate the 'diagonals' covariance
     covg3  =covgammiosl(th1,params,3,ifinv);
     cov3c=varianceofestimates(th1,params,covg3,ifinv);
     try
       % Store the unique values in [s2s2 nunu rhrh s2nu s2rh nurh] order
       covempt(ind,:)=trilos(covemp);
       cov1t(ind,:)  =trilos(cov1c);
       cov2t(ind,:)  =trilos(cov2c);
       cov3t(ind,:)  =trilos(cov3c);
     catch
       keyboard
     end
   end
   cmnylim=[minmax([abs(covempt) abs(cov1t) abs(cov2t) abs(cov3t)])];
   cmnylimrng=round(log10(cmnylim),0);
   for jnd=1:6
     axes(ah(jnd))
     % Plot reference diagonals as light gray lines on bottom for each power of
     % 10 between cmnylim
     hold on
     rnd=0;
     if cs==1
       for ind=cmnylimrng(1):cmnylimrng(2)
         rnd=rnd+1;
         refd(ind)=loglog(dydxs,10^(ind)./NyNxs,'Color',[0.8 0.8 0.8]);
       end
     elseif cs==2
       for ind=cmnylimrng(1):cmnylimrng(2)
         rnd=rnd+1;
         refd(ind)=loglog(dydxs,10^(ind)./gslen,'Color',[0.8 0.8 0.8]);
       end
     end
     % Plot the absolute covariance so that negative correlations do not 
     % mislead
     loglog(dydxs,abs(covempt(:,jnd)),'ko','DisplayName','empirical')
     loglog(dydxs,abs(cov1t(:,jnd))  ,'^','MarkerEdgeColor',clrs{1},...
            'MarkerFaceColor',clrs{1},'DisplayName','sample')
     loglog(dydxs,abs(cov2t(:,jnd))  ,'s','MarkerEdgeColor',clrs{2},...
            'MarkerFaceColor',clrs{2},'DisplayName','dftmtx')
     loglog(dydxs,abs(cov3t(:,jnd))  ,'*','MarkerFaceColor',clrs{3},...
            'MarkerFaceColor',clrs{3},'DisplayName','diagonal')
   end
   for jnd=1:6
     axes(ah(jnd));
     % Fit a slope to the covariance calculations as a function of grid length
     % in log-log
     mbfemp(jnd,:) =polyfit(log10(dydxs),log10(abs(covempt(:,jnd))),1);
     mbfcov1(jnd,:)=polyfit(log10(dydxs),log10(abs(cov1t(:,jnd))),1);
     mbfcov2(jnd,:)=polyfit(log10(dydxs),log10(abs(cov2t(:,jnd))),1);
     mbfcov3(jnd,:)=polyfit(log10(dydxs),log10(abs(cov3t(:,jnd))),1);
     predemp(jnd,:) =polyval(mbfemp(jnd,:), log10(dydxs));
     predcov1(jnd,:)=polyval(mbfcov1(jnd,:),log10(dydxs));
     predcov2(jnd,:)=polyval(mbfcov2(jnd,:),log10(dydxs));
     predcov3(jnd,:)=polyval(mbfcov3(jnd,:),log10(dydxs));
     loglog(dydxs,10.^(predemp(jnd,:)),'Color','k');
     loglog(dydxs,10.^(predcov1(jnd,:)),'Color',clrs{1});
     loglog(dydxs,10.^(predcov2(jnd,:)),'Color',clrs{2});
     loglog(dydxs,10.^(predcov3(jnd,:)),'Color',clrs{3});

     % Annotate the slopes
     % Not enough space for the following:
     % leg1(ind)=legend(sprintf('emp, %s%0.1f}','(\pi\rho)^{',mbfemp(1,1)),...
     %                  sprintf('sam, %s%0.1f}','(\pi\rho)^{',mbfcov1(1,1)),...
     %                  sprintf('dft, %s%0.1f}','(\pi\rho)^{',mbfcov2(1,1)),...
     %                  sprintf('dia, %s%0.1f}','(\pi\rho)^{',mbfcov3(1,1)),...
     %                  'interpreter','tex','BackgroundAlpha',0,'box','off');
     % instead, just report the slope and describe in a figure caption:
     [leg1(jnd),legic]=legend('','','',...
                      sprintf('empirical, %0.1f',abs(mbfemp(jnd,1))),...
                      sprintf('sample, %0.1f',abs(mbfcov1(jnd,1))),...
                      sprintf('dftmtx, %0.1f',abs(mbfcov2(jnd,1))),...
                      sprintf('diagonal, %0.1f',abs(mbfcov3(jnd,1))),...
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
     title(sprintf('%s',labs{jnd}))
     longticks
     ylim([cmnylim(1)*0.5 cmnylim(2)*2.0])
   end      
    
   if cs==1
    xlb='Spacing ($\Delta x)$';
   elseif cs==2
    % TODO is this true?
    xlb='Length ($\pi\rho / \Delta x$)';
   end
   % Plot the bias of the estimates from the empirical study
   axes(ah(7));
   semilogx(dydxs,squeeze(mean(thhats(:,1,:))-th1(1)),'ko')
   % Add the 0-line where the estimate equals the truth
   yline(0)
   xlabel(xlb,'Interpreter','latex')
   ylabel('$\langle\hat{\theta}\rangle-\theta_0$','interpreter','latex');
   ti(7)=title('\sigma^2','FontWeight','bold','Interpreter','tex');
   longticks
   %
   axes(ah(8));
   semilogx(dydxs,squeeze(mean(thhats(:,2,:))-th1(2)),'ko')
   yline(0)
   xlabel(xlb,'Interpreter','latex')
   ti(8)=title('\nu','FontWeight','bold','Interpreter','tex');
   ylim([-0.1,max(mean(thhats(:,2,:))-th1(2))+0.1])
   longticks
   %
   axes(ah(9));
   semilogx(dydxs,squeeze(mean(thhats(:,3,:))-th1(3)),'ko')
   yline(0)
   xlabel(xlb,'Interpreter','latex')
   ti(9)=title('\rho','FontWeight','bold','Interpreter','tex');
   longticks
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
     sprintf('covgammiosl_demo8_paramcovgridsize_s2%inu%irh%i_%s.eps',...
       th1.*[1 10 1],date),'epsc')
end
