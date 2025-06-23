function varargout=mleosl(Hx,thini,params,algo,bounds,aguess,ifinv,xver,cm)
% [thhat,covFHhJ,lpars,scl,thini,params,Hk,k]=...
%          MLEOSL(Hx,thini,params,algo,bounds,aguess,ifinv,xver,cm)
%
% Maximum-likelihood estimation for univariate Gaussian
% multidimensional fields with isotropic Matern covariance
% See Olhede & Simons (2013), doi: 10.1093/gji/ggt056.x
% See Guillaumin et al. (2022), doi: 10.1111/rssb.12539
%
% INPUT:
%
% Hx      Real-valued column vector of unwrapped spatial-domain quantities 
% thini   An unscaled starting guess for the parameter vector with elements:
%         [          s2 nu rho], see SIMULOSL. If you leave this value
%         blank, then you will work from the perturbed "aguess" 
% params  A parameter structure with constants assumed known, see SIMULOSL
%         [dydx NyNx blurs kiso] in the units of 
%         m (2x), "nothing" (3x), rad/m, namely, in order:
%         blurs 0 No wavenumber blurring
%               1 No wavenumber blurring, effectively
%               N Fejer convolutional  BLUROS  on an N-times refined grid
%              -1 Fejer multiplicative BLUROSY using exact procedure
%             Inf Error -> Only for SIMULOSL to use SGP invariant embedding
%         kiso   wavenumber beyond which we are not considering the likelihood
%         quart 1 quadruple, then QUARTER the spatial size
%               0 size as is, watch for periodic correlation behavior
%         taper 0 there is no taper near of far, same as 1
%               1 it's a unit taper, implicitly, same as 0
%               OR an appropriately sized taper with proper values 
%                  (1 is yes and 0 is no and everything in between)
% algo    'unc' uses FMINUNC for unconstrained optimization [default]
%         'con' uses FMINCON with positivity constraints
%         'dsm' uses FMINSEARCH for downhill simplex method, derivative-free
%         'klose' simply closes out a run that got stuck [defaulted when needed]
% bounds   A cell array with those positivity constraints [defaulted]
% aguess   A parameter vector [s2 nu rho] that will be used in
%          simulations for demo purposes, and on which "thini" will be
%          based if that was left blank. If "aguess" is blank, there is
%          a default. If "thini" is set, there is no need for "aguess"
% ifinv    ordered inversion flags for [s2 nu rho], e.g., [1 0 1]:
%          only minimizes the log-likelihood for the 1-flagged parameters, with  
%          the 0-flagged fixed to the value provided in thini, directly or as default;
%          this has the effect of speeding up the estimation procedure but may
%          not be appropriate in the case of real data with various unknowns
% xver     0 Minimal output, no extra verification steps
%          1 Conduct extra verification steps
% cm       1, 2, 3 also launches COVTHOSL with this number as method, and 
%          0 does not (because doing it takes time)
%
% OUTPUT:
%
% thhat    The maximum-likelihood estimate of the vector [scaled]:
%          [s2 nu rho], in units of variance, "nothing", and distance, see SIMULOSL
% covFHhJ  The covariance estimates:
%          covFHhJ{1} from Fisher matrix AT the estimate [FISHIOSL] (eq. 139)
%          covFHhJ{2} from analytical Hessian matrix AT the estimate [HESSIOSL] (eq. 133)
%          covFHhJ{3} from numerical Hessian matrix NEAR the estimate [FMINUNC/FMINCON]
%          covFHhJ{4} from the full formula (eq. 138) at the estimate [COVTHOSL] if you have it
% lpars    The logarithmic likelihood and its derivatives AT or NEAR the estimate
%          lpars{1} the numerical logarithmic likelihood [FMINUNC/FMINCON]
%          lpars{2} the numerical scaled gradient, or score [FMINUNC/FMINCON]
%          lpars{3} the numerical scaled second derivative, or Hessian [FMINUNC/FMINCON]
%          lpars{4} the exit flag of the FMINUNC/FMINCON procedure [bad if 0]
%          lpars{5} the output structure of the FMINUNC/FMINCON procedure
%          lpars{6} the options used by the FMINUNC/FMINCON procedure
%          lpars{7} any bounds used by the  FMINUNC/FMINCON procedure
%          lpars{8} the residual moment statistics used for model testing 
%          lpars{9} the predicted variance of lpars{8}(3) under the null hypothesis
% scl      The scaling that applies to THHAT and THINI
% thini    The starting guess used in the optimization procedure [scaled]
% params   The known constants used inside, see above under INPUT
% Hk       The spectral-domain version of the spatial-domain vector Hx
% k        The wavenumbers on which the estimate is actually based
%
% NOTE: 
%
% A program like EGGERS5 runs 'demo1' in an SPMD loop. Files are opened in
% append mode, except "thzro", which only reflects one lab in that case.
% Writing wires could get cross-checked, messing up the "diagn" files.
%
% EXAMPLE:
%
% % 2015 - convolutional blurring all the way
% p.quart=0; p.blurs=2; p.kiso=NaN; clc; [Hx,th,p]=simulosl([],p,1);
% p.blurs=2; mleosl(Hx,[],p,[],[],[],[],1,1);
%
% % 2018 - exact blurring all the way
% p.quart=0; p.blurs=-1; p.kiso=NaN; clc; [Hx,th,p]=simulosl([],p,1);
% p.blurs=-1; mleosl(Hx,[],p,[],[],[],[],1,0);
%
% % 2025 - circulant embedding and exact blurring
% p.quart=0; p.blurs=Inf; p.kiso=NaN; clc; [Hx,th,p]=simulosl([],p,1);
% p.blurs=-1; mleosl(Hx,[],p,[],[],[],[],1,1);
%
% You can stick in partial structures, e.g. only specifying params.kiso
%
% Fix any of the parameters to aguess by setting ordered value of ifinv to 0 
% mleosl(Hx,[],p,[],[],[],[1 0 0],1)
%
% Invert for the squared exponential case for nu fixed to Inf 
% Hx=simulosl([1e6 Inf 1e4],p,1);
% mleosl(Hx,[],p,[],[],[1e6 Inf 1e4],[1 0 1],1)
%
% Perform a series of N simulations centered on th0 with different p's
% mleosl('demo1',N,th,p)
%
% Statistical study of a series of simulations using MLEPLOS
% mleosl('demo2','18-May-2025')
%
% Covariance study of a series of simulations using COVPLOS
% mleosl('demo4','14-Oct-2023')
%
% One simulation and a chi-squared plot using MLECHIPLOS
% mleosl('demo5',th,p) % This should be as good as
% blurosy('demo2',p.NyNx,[],[],1) % for the same th
%
% Estimation of one realization for 3-parameter and 2-parameter optimization
% mleosl('demo6')
%
% Demo series for a fixed parameter
% mleosl('demo1',N,th,p,[],[],th,[1 0 1]);
% mleosl('demo2',date,[],[],[],[],[1 0 1]);
% mleosl('demo4',date,[],[],[],[],[1 0 1]);
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
%
% Last modified by fjsimons-at-alum.mit.edu, 05/30/2025
% Last modified by olwalbert-at-princeton.edu, 05/30/2025

if ~isstr(Hx)
  defval('algo','unc')
  % The necessary strings for formatting, see OSDISP and OSANSW
  str0='%18s';
  str1='%13.0e '; 
  str2='%13.0f %13.2f %13.0f';
  str2='%13.3g %13.3g %13.5g';
  str3s='%13s ';

  % Supply the needed parameters, keep the givens, extract to variables
  fields={               'dydx','NyNx','blurs','kiso','quart','taper'};
  defstruct('params',fields,...
	    {                      [20 20]*1e3,sqrt(length(Hx))*[1 1],-1,NaN,0,0});
  struct2var(params)

  % You cannot call MLEOSL with params.blurs=Inf, since that's for
  % SIMULOSL only, we reset for the inversion only inside LOGLIOSL

  % These bounds are physically motivated...
  if strcmp(algo,'con') || strcmp(algo,'dsm')
    % Parameters for FMINCON in case that's what's being used, which is recommended
    defval('bounds',{[],[],... % Linear inequalities
                     [],[],... % Linear equalities
                     [0.1 0.15 sqrt(prod(dydx))],... % Lower bounds
                     [100 8.00 max(2.5e5,min(dydx.*NyNx))],... % Upper bounds
                     []}); % Nonlinear (in)equalities
  else
    bounds=[];
  end

  % Being extra careful or not?
  defval('xver',1)

  % Invert for [s2 nu rh]?
  defval('ifinv',[1 1 1])

  % The parameters used in the simulation for demos, or upon which to base "thini"
  % Check Vanmarcke 1st edition for suggestions on initial rho, very important

  % An alternative guessing approach                                            
  altguess=0;                                                                   
  if altguess                                                                   
      disp('in development: altguess=1; flag off when done')
      % Approximate variance parameter (e.g., eq 32 of Methodology paper)       
      s2guess=nanvar(Hx);                                                       
      % We can estimate the range parameter from the analytical definition of 
      % the spectral density (Eq. 72 of S&O 2013), the empirical periodogram at
      % the origin, and our guess for s2
      [k,~,~,kxv,kyv]=knums(params);%,1);

      % For an idealized situation where we know th0 (e.g., used for testing)
      % th0=[1e6 1.5 2e4];
      % Sk=v2s(maternos(k(:),th0));
      % prdgrm=Sk;

      % For the empirical periodogram (consider the taper or not?)
      % Hk2=v2s(abs(tospec(Tx(:).*Hx(:,1),params)/(2*pi)).^2,params);
      Hk2=v2s(abs(tospec(Hx(:,1),params)/(2*pi)).^2,params);
      prdgrm=Hk2;

      % At the origin
      S0=prdgrm(~k(:));                                                         

      % Assuming a 2-dimensional field     
      rhoguess=2*sqrt(S0/(s2guess*pi));                                        
                                                                                
      % Selection of initial nu value inspired by Sykulski et al 2019 Biometrika
      % (doi:10.1093/biomet/asy071) where initialization values for spectral 
      % parameters of a correlation structure are determined via ordinary least
      % squares for 1-D data. I attempted to find slope of best-fit plane for a 
      % quadrant of the 2-D periodogram with OLS and SVD, but the results were 
      % poor. Instead, finding the slope of the radial average allows for an OK 
      % guess, with the exception of a dependency on rho that would not have
      % have been a challenge in the Biometrika paper based on their analytical 
      % parameterization. For now, scaling nuguess by the scaling predicted for
      % the range seems sufficient enough to get the inversion ticking, which is
      % all we are after anyway.
      rhscl=max(1,10.^round(log10(abs(rhoguess))));
      % Assemble a vector of low wavenumbers where the slope of the radial
      % average is near-constant; the current range was selected by eye from 
      % the commented out plotting routine below
      prdgrm(k>max(abs(kxv)))=NaN;
      [~,rav]=radavg(prdgrm);
      [~,kav]=radavg(k);
      K=length(kav);                                                             
      rkmin=max(4,floor(K*0.05));
      rkmax=2*rkmin;                                      
      ktrm=kav(rkmin:rkmax);
      rtrm=rav(rkmin:rkmax);  

      % Least squares for slope parameter
      G = [ones(prod(size(rtrm)),1) ktrm];
      b=inv(G'*G)*G'*rtrm(:); 
      bscl=10.^(round(mean(log10(rtrm))));
      nuguess=max(0.5,-1/2*real(b(2)/2/bscl/rhscl));

      % LS with log of perturbed data; doesn't perform as well
      % Glog = [ones(prod(size(rtrm)),1) log(ktrm)];
      % blog=inv(Glog'*Glog)*Glog'*log(rtrm(:)); 
      % nuguesslog=max(0.5,-1/2*real(blog(2)/2)) 

      % Plot the periodogram and the traces we select to see how well our
      % radial average and its linear fit approximate the periodogram
      % figure();imagesc(prdgrm);
      % figure();plot(kav,rav);hold on;plot(ktrm,rtrm)
      % plot([kav(rkmin) kav(rkmin)],[1e-20 max(rtrm)],...
      %      [kav(rkmax) kav(rkmax)],[1e-20 max(rtrm)]);
      % plot(ktrm,b(1)+b(2)*ktrm)
      defval('aguess',[s2guess nuguess rhoguess]);
  else
      defval('aguess',[nanvar(Hx) 2.0 sqrt(prod(dydx.*NyNx))/pi/2/10]);
  end

  % Scale the parameters by this factor; fix it unless "thini" is supplied
  defval('scl',10.^round(log10(abs(aguess))));

  % Unless you supply an initial value, construct one from "aguess" by perturbation
  nperturb=0.25;
  % So not all the initialization points are the same!!
  defval('thini',abs((1+nperturb*randn(size(aguess)).*ifinv).*aguess))
  % If you brought in your own initial guess, you need an appropriate new scale
  if ~isempty(inputname(2)) || any(aguess~=thini)
    scl=10.^round(log10(abs(thini)));
    disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str1,size(scl))),...
		 'Scaling',scl))
  end
  disp(sprintf(sprintf('%s : %s',str0,str2),...
	       'Starting theta',thini))
  if strcmp(algo,'con')
    disp(sprintf(sprintf('\n%s : %s',str0,str2),...
                 'Lower bounds',bounds{5}))
    disp(sprintf(sprintf('%s : %s',str0,str2),...
                 'Upper bounds',bounds{6}))
  end

  % Now apply the scale so the minimization doesn't get into trouble; for the 
  % special case of nu->Inf, we must enforce that scl(2)==1
  if thini(end-1)==Inf; scl(2)=1; end
  thini=thini./scl;
 
  % Analysis taper
  if length(taper)==1 && (taper==0 || taper==1)
      Tx=1;
  else
      % Ones and zeros as suitable for BLUROSY
      Tx=taper;
  end

  % Create the appropriate wavenumber axis
  k=knums(params);
  
  % We could get into the habit of never involving the zero-wavenumber
  knz=(~~k);

  % Always scale the data sets but don't forget to reapply at the very end
  shat=nanstd(Hx(:,1));
  % Scale the data; don't reorder the next three lines!
  Hx(:,1)=Hx(:,1)./shat;
  % Rescale the initial value so the output applies to both THHAT and THINI
  thini(1)=thini(1).*scl(1)/shat^2;
  % And with these new scalings you have no more business for the first scale
  % though for the derived quantity you need to retain them
  sclh=[shat^2 scl(2:3)];
  matscl=[sclh(:)*sclh(:)'];
  % Every next occurence has theta(1) compared to scaled Hk so that makes scl(1)=1
  scl(1)=1;
  matscl1=[scl(:)*scl(:)'];

  % If we have the special case of nu->Inf, we will make a special matrix scale
  if isinf(thini(end-1)) && ifinv(end-1)==0
    matsclInf=matslice(matscl,[1 0 1]);
  end

  % Always demean the data sets - think about deplaning as well?
  Hx(:,1)=Hx(:,1)-nanmean(Hx(:,1));

  % Turn the tapered observation vector to the spectral domain
  % Watch the 2pi in SIMULOSL
  Hk(:,1)=tospec(Tx(:).*Hx(:,1),params)/(2*pi);

  % Account for the size here? Like in SIMULOSL and checked in BLUROSY
  % See BLUROSY and how to normalize there, maybe take values of Tx?
  if size(Tx)~=1
      % This adjusts for the size of an explicit taper
      Hk(:,1)=Hk(:,1)/sqrt(sum(Tx(:).^2))*sqrt(prod(params.NyNx));
  end
  
  NN=200;
  % And now get going with the likelihood using Hk
  % [ off|iter|iter-detailed|notify|notify-detailed|final|final-detailed ] 
  % Should probably make the tolerances relative to the number of k
  % points? Or watch at least how these gradients size to the scaled lik 
  options=optimset('GradObj','off','Display','off',...
		   'TolFun',1e-11,'TolX',1e-11,'MaxIter',NN,...
		   'LargeScale','off');
  % The 'LargeScale' option goes straight for the line search when the
  % gradient is NOT being supplied.

  % Set the parallel option to (never) use it for the actual optimization
  % Doesn't seem to do much when we supply our own gradient
  % options.UseParallel='always';

  % The number of parameters that are being solved for
  np=length(thini);
  % The number of unique entries in an np*np symmetric matrix
  npp=np*(np+1)/2;

  if xver==1 && blurs>-1 && blurs<2 
    % Using the analytical gradient in the optimization is not generally a good
    % idea but if the likelihoods aren't blurred, you can set this option to
    % 'on' and then let MATLAB verify that the numerical calculations match
    % the analytics. According to the manual, "solvers check the match at a
    % point that is a small random perturbation of the initial point". My
    % own "disp" output (further below) provides comparisons at the estimate
    % even when the option below is set to "off" and we don't use it for any
    % aspect of the optimization. If you should also try this for blurred
    % systems (remove part of the condition above), you will fail the
    % test and the whole thing will come to a halt. So after doing this
    % interactively a few times, I've been setting the below to "off". 
    options.GradObj='off';
    % Leave the below "on" since it's inconsequential when the above is "off"
    options.DerivativeCheck='on';
  end

  % And find the MLE! Work on scaled parameters
  try
    switch algo
     case 'unc'
      % disp('Using FMINUNC for unconstrained optimization of LOGLIOSL')
       if any(ifinv~=[1 1 1])
           % We want to optimize for only the parameters requested, taking the
           % value of thini for fixed parameters
           % Put the non-inverted-for parameters up front
           [~,tmpidx]=sort(ifinv);
           % But remember how to put them back in order for later
           [~,invidx]=sort(tmpidx);

           t0=clock;
           % It's the FMINUNC that benefits from the scaling, but LOGLIOSL is in units
           [thhat,logli,eflag,oput,~,~]=...
               fminunc(@(theta) logliosl(k,scl.*indeks([thini(~ifinv) theta],invidx),...
                                         params,Hk,xver),...
                       thini(~~ifinv),options);
           ts=etime(clock,t0);

           % The estimate of theta that we report should include the optimized
           % and fixed parameters
           thhat=indeks([thini(~ifinv) thhat],invidx);
           % Calculate the numerical gradient and Hessian given all three
           % parameters
           derivopts=optimset('MaxIter',0,'MaxFunEvals',0,'Display','off');
           [~,~,~,~,grd,hes]=...
               fminunc(@(theta) logliosl(k,scl.*theta,...
                                         params,Hk,0),...
                       thhat,derivopts);
           % In the special case of nu->Inf, we will retain a slice of grd and 
           % hes to avoid singularities when we later calculate the numerical 
           % covariance approximations
           if isinf(thini(end-1)) && ifinv(end-1)==0
             grdInf=matslice(grd,ifinv);
             hesInf=matslice(hes,ifinv);
           end
       else
           t0=clock;
           [thhat,logli,eflag,oput,grd,hes]=...
	       fminunc(@(theta) logliosl(k,scl.*theta,...
                                         params,Hk,xver),...
		       thini,options);
           ts=etime(clock,t0);
       end
     case 'con'
      % New for FMINCON
      options.Algorithm='active-set';
      % disp('Using FMINCON for constrained optimization of LOGLIOSL')
      t0=clock;
      % See M. K. Stein p. 173 when differentiability parameter maxes out
      % Also check for when this crucial parameter hits the lower bound.
      % Hitting the bounds for a parameter is relative, say within 10%.
      thhat(2)=bounds{6}(2); nwh=4; nwi=0; hitit=thhat(2)/10;
      while [bounds{6}(2)-thhat(2)<hitit ...
             || thhat(2)-bounds{5}(2)<hitit] ...
             && nwi<nwh
	  nwi=nwi+1;
	  thisthini=thini;
          if nwi>1
              disp(sprintf(...
                  '\nHit the wall on differentiability... trying again %i/%i\n',...
                  nwi,nwh))
          end
          if any(ifinv~=[1 1 1])
              % Only optimize for the parameters requested
              for bdx=1:length(bounds)
                % For the lower and upper bounds of the parameter set
                  if ~isempty(bounds{bdx})
                      % Only scale the parameter bounds actually inverted for
                    subbounds{bdx}=bounds{bdx}(logical(ifinv))./scl(logical(ifinv));
                else
                    % This remains empty
                    subbounds{bdx}=[];
                  end
              end

              % Put the non-inverted-for parameters up front
              [~,tmpidx]=sort(ifinv);
              % But remember how to put them back in order for later
              [~,invidx]=sort(tmpidx);

              [thhat,logli,eflag,oput,lmd,~,~]=...
                  fmincon(@(theta) logliosl(k,scl.*indeks([thini(~ifinv) theta],invidx),...
                                            params,Hk,xver),...
                          thini(~~ifinv),...
                          subbounds{1},subbounds{2},subbounds{3},subbounds{4},...
                          subbounds{5},subbounds{6},subbounds{7},...
                          options);
              % The estimate of theta should include all three parameters, whether
              % we fixed or optimized for them
              thhat=indeks([thini(~ifinv) thhat],invidx);
              % Whatever the incoming bounds were, make the lower and upper bound of
              % the non-inverted-for parameters equal to their fixed value
              bounds{5}=indeks([thini(~ifinv)./scl(~ifinv) subbounds{5}],invidx);
              bounds{6}=indeks([thini(~ifinv)./scl(~ifinv) subbounds{6}],invidx);

              % Calculate the numerical gradient and Hessian from all three
              % parameters; keep in mind that FMINCON estimates numerical
              % Hessian at the next-to-last iteration from the Lagrangian 
              % rather than the objective function directly; may present
              % inaccuracies
              derivopts=optimset('MaxIter',0,'MaxFunEvals',0,'Display','off');
              [~,~,~,~,~,grd,hes]=...
                  fmincon(@(theta) logliosl(k,scl.*theta,...
                                            params,Hk,0),...
                          thhat,...
                          bounds{1},bounds{2},bounds{3},bounds{4},...
                          bounds{5},bounds{6},bounds{7},...
                          derivopts);
              % In the special case of nu->Inf, we will retain a slice of grd and 
              % hes to avoid singularities when we later calculate the numerical 
              % covariance approximations
              if isinf(thini(end-1)) && ifinv(end-1)==0
                  grdInf=matslice(grd,ifinv);
                  hesInf=matslice(hes,ifinv);
              end
          else
              [thhat,logli,eflag,oput,lmd,grd,hes]=...
                  fmincon(@(theta) logliosl(k,scl.*theta,...
                                            params,Hk,xver),...
                          thini,...
                          bounds{1},bounds{2},bounds{3},bounds{4},...
                          bounds{5}./scl,bounds{6}./scl,bounds{7},...
                          options);
              % Try resetting the offending parameter nu by a serious kick
              thini(end-1)=thini(end-1)/[1+1/4-rand/2];
          end
          % And the others, switching the relationship between sigma^2 and rho
          thini(1)=thini(1)*rand;
          thini(3)=thini(3)/rand;
      end
      % You've left the loop, so you've used the last thini
      thini=thisthini;
      if nwi==nwh
        % You haven't been able to do it within the bounds for nu, relax
        % the bounds or flag the result, or trim it at the end... error!
        warning('Solution hugs the bound for NU, perhaps uncomfortably')
      end
      ts=etime(clock,t0);
     case 'dsm'
         % FMINSEARCH minimizes LOGLIOSL objective function according to the 
         % Nelder-Mead simplex algorithm by first generating a simplex with n+1
         % vertices about the n-dimensional thini, and then iteratively 
         % modifying the simplex through reflection, expansion, and contraction
         % of function value-ordered vertices; must satisfy both the function
         % tolerance and step tolerance stopping criteria to exit. This is a 
         % gradient-free solver so convergence to local minimum not guaranteed. 
         % FMINSEARCH was used by the Biometrika 2019 publication 
         % doi:10.1093/biomet/asy071
         if any(ifinv~=[1 1 1])
           % We want to optimize for only the parameters requested, taking the
           % value of thini for fixed parameters
           % Put the non-inverted-for parameters up front
           [~,tmpidx]=sort(ifinv);
           % But remember how to put them back in order for later
           [~,invidx]=sort(tmpidx);

           t0=clock;
           options=[];
           [thhat,logli,eflag,oput]=...
             fminsearch(@(theta) logliosl(k,scl.*indeks([thini(~ifinv) theta],invidx),...
                                          params,Hk,xver),...
                        thini(~~ifinv),options);
           ts=etime(clock,t0);
           % The estimate of theta should include all three parameters, whether
           % we fixed or optimized for them
           thhat=indeks([thini(~ifinv) thhat],invidx);
        else
           t0=clock;
           options=[];
           [thhat,logli,eflag,oput]=...
	       fminsearch(@(theta) logliosl(k,scl.*theta,...
                                            params,Hk,xver),...
		          thini,options);
           ts=etime(clock,t0);
        end
        % While FMINSEARCH is gradient-free, we might still want to know
        % what these values are at the estimate. Let's find the numerical
        % gradient and Hessian at the parameter estimate following our
        % approach for doing so in the fixed parameter case for FMINUNC
        derivopts=optimset('MaxIter',0,'MaxFunEvals',0,'Display','off');
        
        [~,~,~,~,grd,hes]=...
            fminunc(@(theta) logliosl(k,scl.*theta,...
                                      params,Hk,0),...
                    thhat,derivopts);
        % In the special case of nu->Inf, we will retain a slice of grd and 
        % hes to avoid singularities when we later calculate the numerical 
        % covariance approximations
        if isinf(thini(end-1)) && ifinv(end-1)==0 
           grdInf=matslice(grd,ifinv);
           hesInf=matslice(hes,ifinv);
        end
      case 'klose'
       % Simply a "closing" run to return the options
       lpars{6}=options;
       lpars{7}=bounds;
       % Simply a "closing" run to return the options
       varargout=cellnan(nargout,1,1);
       varargout{end}=lpars;
       return
    end
  catch
    % If something went wrong, exit gracefully
    if ~exist('thhat') || isnan(thhat); thhat=deal(nan(1,3)); end
    varns={thhat,[],[],sclh,thini,params,Hk,k};
    varargout=varns(1:nargout);
    return
  end

  % It is not impossible that a solution is reached which yields a
  % negative rho - which only appears in the square in MATERNOS. But if
  % we're going to calculate (approximately blurred) analytical
  % gradients and Hessians (even using exact blurring of the spectral
  % densities), we are going to be using MATERNOSY, which will complain...
  if thhat(1)<0
    error(sprintf('%s Negative variance',upper(mfilename)))
  end
  if thhat(2)<0
    error(sprintf('%s Negative smoothness',upper(mfilename)))
  end
  if thhat(3)<0
    error(sprintf('%s Negative range',upper(mfilename)))
  end

  % Degrees of freedom for full-wavenumber domain (redundant for real data)
  % Not including the zero wavenumber, since LOGLIOS doesn't either
  df=length(k(~~k))/2;

  % Watch out for singularity or scaling warnings, they are prone to pop up

  % Covariance from FMINUNC/FMINCON's numerical scaled Hessian AT/NEAR estimate,
  % taking special care for the special case of nu->Inf
  if isinf(thini(end-1)) && ifinv(end-1)==0
      covh=inv(hesInf./matsclInf)/df;
  else
      covh=inv(hes./matscl)/df;
  end
  
  if xver==1 & verLessThan('matlab','8.4.0')
    % Try variable-precision arithmetic?
    vh=sym('vh',[np np]);
    for index=1:prod(size(vh))
      vh(index)=sprintf('%0.16e/%0.1e ',hes(index),matscl(index));
    end
    % Could try even more digits with VPA but in the end it didn't all seem
    % to matter much
    vcovh=inv(vh)/df;
  end
    
  % Fisher matrix AT the estimate, and covariance derived from it
  [F,covF]=fishiosl(k,sclh.*thhat,params,xver);
  if isinf(thini(end-1)) && ifinv(end-1)==0
    FInf=matslice(F,ifinv);
    covF=inv(FInf)/df;
  end

  % Analytic (poorly blurred) Hessian AT the estimate, and derived covariance
  [H,covH]=hessiosl(k,scl.*thhat,params,Hk,xver);
  scls=[shat^2 1 1];
  matscls=[scls(:)*scls(:)'];
  covH=inv(-H./matscls)/df;

  % For the special case of nu->Inf, we will neglect terms involving nu
  if isinf(thini(end-1)) && ifinv(end-1)==0
      HInf=matslice(H,ifinv);
      covH=inv(-HInf)/df;
  end

  % FJS how about a step further, use F-1 H F-T to get any influence at all
  % Does Arthur use the average variance of the gradient here somewhere
  if isinf(thini(end-1)) && ifinv(end-1)==0
    covFHF=inv(FInf)*[-HInf]*inv(FInf)/df;
    covFhF=inv(FInf)*[hesInf./matsclInf]*inv(FInf)/df;
  else
    covFHF=inv(F)*[-H./matscls]*inv(F)/df;
    covFhF=inv(F)*[hes./matscl]*inv(F)/df;
  end

  % Calculate the variance of the estimates using COVGAMMIOSL via 
  % (1) gradient sampling (always works, robust, fast)
  % (2) full dftmtx (fast, but very memory intensive)
  % (3) diagonals method (slow, ultimately)
  % (0) don't even do it
  defval('cm',0);
  covFJF=covthosl(sclh.*thhat,params,cm,ifinv);
  % This should ensure that COVFJF is a 3x3 whether we inverted for all
  % parameters or not; only needed for display purposes later
  covFJF=matslice(covFJF,ifinv,-1);

  % Analytical calculations of the gradient and the Hessian poorly represent
  % the blurring (though it's much better than not trying at all), and thus,
  % are expected to be close to numerical results only without blurring
  if xver==1
    % Analytic (poorly blurred) gradient, scaled for numerical comparison
    gros=gammiosl(k,scl.*thhat,params,Hk,xver).*scl(:);

    % Compare the analytic Hessian with the numerical Hessian and with
    % the Hessian expectation, which is the Fisher, at the estimate, and
    % compare the analytic gradient with the numerical gradient
    str3=repmat('%13g ',1,npp);
    str4=repmat('%13g ',1,np);
    disp(sprintf('%s',repmat('_',119,1)))
    disp(sprintf('\n%16s\n','At the ESTIMATE:'));
    disp(sprintf(sprintf('    Log-likelihood : %s',str3),logli))

    disp(sprintf(...
	['\nThe numerical derivatives are usually at the penultimate iteration:']))
    if params.blurs~=0
      disp(sprintf(...
          ['\nWith blurring, the comparisons below are necessarily inexact:']))
    end
    disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3s,1,np)),...
	       ' ','ds2','dnu','drho'))
    disp(sprintf(sprintf(' Numericl Gradient : %s',str4),grd))
    disp(sprintf(sprintf(' Analytic Gradient : %s',str4),gros))

    disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3s,1,npp)),...
	       ' ','(ds2)^2','(dnu)^2','(drho)^2','ds2dnu','ds2drho','dnudrho'))
    disp(sprintf(sprintf(' Numerical hessian : %s',str3),trilos(hes)))
    disp(sprintf(sprintf(' Analyticl Hessian : %s',str3),trilos(-H.*matscl1)))
    disp(sprintf(sprintf(' Analytical Fisher : %s',str3),trilos( F.*matscl)))

    disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3s,1,npp)),...
	       ' ','C(s2,s2)','C(nu,nu)','C(rho,rho)','C(s2,nu)','C(s2,rho)','C(nu,rho)'))
    disp(sprintf(sprintf(' Cov (Numer hess.) : %s',str3),trilos(covh)))
    disp(sprintf(sprintf(' Cov (Analy Hess.) : %s',str3),trilos(covH)))
    disp(sprintf(sprintf(' Cov (Analy Fish.) : %s',str3),trilos(covF)))
    disp(sprintf(sprintf(' Cov ( FishhFish.) : %s',str3),trilos(covFhF)))
    disp(sprintf(sprintf(' Cov ( FishHFish.) : %s',str3),trilos(covFHF)))
    if sum(covFJF(:))~=0
        disp(sprintf(sprintf(' Cov ( FishJFish.) : %s',str3),trilos(covFJF)))
    else
        disp(sprintf(sprintf(' Cov ( FishJFish.) : %s',str3),zeros(npp,1)))
    end

    % disp(sprintf('%s\n',repmat('_',119,1)))
    disp(' ')
    disp(sprintf(sprintf('%s : %s ',str0,str2),...
	         'Numer Hessi std',sqrt(diag(covh))))
    disp(sprintf(sprintf('%s : %s ',str0,str2),...
	         'Analy Hessi std',sqrt(diag(covH))))
    disp(sprintf(sprintf('%s : %s\n ',str0,str2),...
	         'Anal Fisher std',sqrt(diag(covF))))
    disp(sprintf(sprintf('%s : %s ',str0,str2),...
	         ' FishhFish. std',sqrt(diag(covFhF))))
    disp(sprintf(sprintf('%s : %s ',str0,str2),...
	         ' FishHFish. std',sqrt(diag(covFHF))))
    disp(sprintf(sprintf('%s : %s\n ',str0,str2),...
	         ' FishJFish. std',sqrt(diag(covFJF))))
  end

  % Talk!
  disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3s,1,np)),...
	       ' ','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str2),...
	       'Estimated theta',sclh.*thhat))
  disp(' ')
  if xver==0 || xver==1
    disp(sprintf('%8.1fs per %i iterations or %5.1fs per %i function counts',...
                 ts/oput.iterations*100,100,ts/oput.funcCount*1000,1000))
    % disp(sprintf('%s\n',repmat('_',119,1)))
  end

  % Here we compute the moment parameters and recheck the likelihood
  [L,~,Hagain,momx,vr]=logliosl(k,scl.*thhat,...
                                params,Hk,xver);
  diferm(L,logli)
  try
      diferm(Hagain,H)
  end

  % Reorganize the output into cell arrays
  covFHhJ{1}=covF;
  covFHhJ{2}=covH;
  covFHhJ{3}=covh;
  covFHhJ{4}=covFJF;
  % Likelihood attributes
  lpars{1}=logli;
  lpars{2}=grd;
  lpars{3}=hes;
  lpars{4}=eflag;
  lpars{5}=oput;
  lpars{6}=options;
  lpars{7}=bounds;
  lpars{8}=momx;
  lpars{9}=vr;
  % Generate output as needed
  varns={thhat,covFHhJ,lpars,sclh,thini,params,Hk,k};
  varargout=varns(1:nargout);
elseif strcmp(Hx,'demo1')
  more off
  % Runs a series of simulations. See 'demo2' to display them.
  % If you run this again on the same date, the files THINI and
  % THHAT get appended, but a blank THZERO is created. 
  defval('thini',[]);
  % How many simulations? The SECOND argument, after the demo id.
  N=thini; clear thini
  % What th-parameter set? The THIRD argument, after the demo id
  defval('params',[])
  % If there is no preference, then that's OK, it gets taken care of
  th0=params; clear params
  % What fixed-parameter set? The FOURTH argument, after the demo id
  defval('algo',[])
  % If there is no preference, then that's OK, it gets taken care of
  params=algo; clear algo
  % What algorithm? The FIFTH argument, after the demo id
  defval('bounds',[])
  % If there is no preference, then that's OK, it gets taken care of
  algo=bounds; clear bounds
  % What datum? The SIXTH argument, after the demo id
  defval('aguess',[])
  % If there is no preference, then that's OK, it gets taken care of
  datum=aguess; clear aguess
  % Where to initialize (approx.)? The SEVENTH argument, after the demo id
  defval('ifinv',[])
  % If there is no preference, then that's OK, it gets taken care of
  aguess=ifinv; clear ifinv
  % Which parameters to invert for? The EIGHTH argument, after the demo id
  defval('xver',[1 1 1])
  % If there is no preference, then that's OK, it gets taken care of
  ifinv=xver; clear xver

  % You can't stick in a NINTH argument so you'll have to default 
  defval('xver',0)

  % What you make of all of that if there hasn't been a number specified
  defval('N',500)

  % The number of parameters to solve for
  np=3;
  
  % Open 'thzro', 'thini', 'thhat' and 'diagn' files and return format strings
  [fids,fmts,fmti]=osopen(np,datum);

  % Do it!
  good=0; 
  % Initialize the average Hessian that will be saved by OSWZEROE 
  avhsz=zeros(np,np);

  % Set N to zero to simply close THZERO out
  for index=1:N
    % Simulate data from the same lithosphere, watch the blurring
    [Hx,th0,p,k,Hk]=simulosl(th0,params,xver);

    % Check the dimensions of space and spectrum are right
    difer(length(Hx)-length(k(:)),[],[],NaN)

    % If we are working within either the squared exponential case or a fixed nu
    % case, we will use the provided 'aguess' within our batch calculuations;
    % otherwise, we will take the provided 'th0' as 'aguess'. 
    if isinf(th0(2)); th1=aguess; elseif ifinv==[1 0 1]; th1=aguess; else; th1=th0; end
    % Form the maximum-likelihood estimate, pass on the params, use th0
    % as the basis for the perturbed initial values. Remember hes is scaled.
    t0=clock;
    [thhat,covFHhJ,lpars,scl,thini,p,Hk,k]=mleosl(Hx,[],p,algo,[],th1,ifinv,xver);
    ts=etime(clock,t0);

    % Initialize the THZRO file... note that the bounds may change
    % between simulations, and only one gets recorded here
    if ~any(isnan(thhat)) && index==1 && labindex==1
      oswzerob(fids(1),th0,p,lpars,fmts)
    end

    % If a model was found, keep the results, if not, they're all NaNs
    % Ignore the fact that it may be at the maximum number of iterations
    % e=1

    % IF NUMBER OF FUNCTION ITERATIONS IS TOO LOW DEFINITELY BAD
    itmin=0;
    % A measure of first-order optimality (which in the unconstrained case is
    % the infinity norm of the gradient at the solution). Maybe what it
    % means to be 'good' should be in function of the data size as more
    % precision will be needed to navigate things with smaller variance! At
    % any rate, you want this not too low.
    optmin=Inf;
    % Maybe just print it and decide later? No longer e>0 as a condition.
    % e in many times is 0 even though the solution was clearly found, in
    % other words, this condition IS a way of stopping with the solution
    try
      % Maybe I'm too restrictive in throwing these out? Maybe the
      % Hessian can be slightly imaginary and I could still find thhat
      if isreal([lpars{1} lpars{2}']) ...
	    && all(thhat>0) ...
	    && all(~isnan(thhat)) ...
	    && lpars{5}.iterations > itmin ...
	    && lpars{5}.firstorderopt < optmin
	good=good+1;
	% Build the AVERAGE of the Hessians for printout by OSWZEROE later
	avhsz=avhsz+lpars{3}./[scl(:)*scl(:)'];
	% Reapply the scalings before writing it out
	fprintf(fids(2),fmts{1},thhat.*scl);
	fprintf(fids(3),fmts{1},thini.*scl);
        % We don't compare the second and third outputs of LOGLIOSL since these are
        % analytical, poorly approximately blurred, derivatives, and we be
        % writing the numerical versions. Be aware that covFHh{3} is the
        % current favorite covariance estimate on the parameters!
	% Print optimization results and diagnostics to different file with OSWDIAG
	oswdiag(fids(4),fmts,lpars,thhat,thini,scl,ts,var(Hx),covFHhJ{3})
      end
    end
  end
  % If there was any success at all, finalize the THZRO file 
  % If for some reason this didn't end well, do an N==0 run.

  % Initialize if all you want is to close the file
  if N==0 
    [Hx,th0,p,k]=simulosl(th0,params); 
    good=1; avhsz=avhsz+1; 
    [~,~,lpars]=mleosl(Hx,[],[],'klose');
    oswzerob(fids(1),th0,p,lpars,fmts)
  end
  
  if good>=1 
    % This is the new scaling based on the truth which we use here 
    sclth0=10.^round(log10(th0));

    % This is the AVERAGE of the numerical Hessians, should be closer to the Fisher
    avhsz=avhsz.*[sclth0(:)*sclth0(:)']/good;

    % If you are working within the squared exponential case, or if you ended on
    % a nonsensical estimate, we will have to intervene and forgo the
    % calculation of the Fisher-derived covariance at the truth
    if ~any(isnan(k(:))) && ~isinf(th0(end-1))
        % Now compute the Fisher and Fisher-derived covariance at the truth
        [F0,covF0]=fishiosl(k,th0,params,xver);
        matscl=[sclth0(:)*sclth0(:)'];
    else
        [F0,covF0,matscl,avhsz]=deal(nan(3,3));
    end

    % Of course when we don't have the truth we'll build the covariance
    % from the single estimate that we have just obtained. This
    % covariance would then be the only thing we'd have to save.
    if labindex==1
      oswzeroe(fids(1),sclth0,avhsz,good,F0.*matscl,covF0,fmti)
    end
  end
  
  % Put both of these also into the thzro file 
  fclose('all');
elseif strcmp(Hx,'demo2')
  defval('thini',[]);
  datum=thini;
  defval('datum',date)
  defval('ifinv',[1 1 1]);

  % The number of parameters to solve for
  np=3;

  % Looks like more trimming is needed for 'con' rather than 'unc'
  trims=100;
  % Load everything you know about this simulation
  [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(datum,trims);

  defval('xver',0)

  if xver==1
      % Report the findings of all of the moment parameters
      disp(sprintf('\nm(m(Xk)) %f m(v(Xk)) %f\nm(magic) %f v(magic) %f',...
	           mean(momx),var(momx(:,end))))
  end
  
  % Plot it all
  figure(1)
  fig2print(gcf,'landscape')
  clf

  % We feed it various things and it calculates a bunch more
  [ah,ha,yl,xl,tl,ps,ms,o1,o2,ep,ec,axlm,covth,ttt]=...
      mleplos(thhats,th0,covF0,covavhs,covXpix,[],[],p,...
                  sprintf('MLEOSL-%s',datum),thpix,ifinv);
  
  % Return some output, how about just the empirically observed
  % means and covariance matrix of the estimates, and the number of
  % reported trials
  nobs=size(thhats,1);
  mobs=mean(thhats,1);
  cobs=cov(thhats);
  varns={cobs,mobs,nobs,th0,p,momx};
  varargout=varns(1:nargout);

  % Print the figure!
  disp(' ')
  % Figure not printing consistently; force MATLAB to revisit the figures drawn
  % prior to save
  figure(2);  figna=figdisp([],sprintf('%s_%s_2',Hx,datum),[],1);
  disp('pause to save')
  pause(3)
  figure(1); figna=figdisp([],sprintf('%s_%s_1',Hx,datum),[],1);

  % Now that the original figure 2 has been saved, modify the figure by removing
  % the points and the cross, and by adding the likelihood from MLELCONTOSL

  flag=3;
  if flag==3
      % Calculate the loglihood contours for an estimate within 
      % MLELCONTOSL -- we may specify the range of parameter values as an input, 
      % but discretization of the grid must be manually modified within 
      % MLELCONTOSL as 'fine'
      % Set the parameter ranges based on the existing axes limits with scaling
      % consistent with nu and rho
      thcont=[axlm(1,:);...
              axlm(2,:).*10^round(log10(th0(2)));...
              axlm(3,:).*10^round(log10(th0(3)))];

      % Set the covariance matrix from which the number of standard deviations
      % out from the mean observation is calculated for contour placement
      covcont=cobs; %covth;

      % Run in parallel? Number of workers assessed as one less number of
      % physical cores.
      runinpar=1;

      % Calculate the grid and objects used for contouring with MLELCONTOSL, all 
      % organized by cross-parameter pairings in the usual order: 
      % s2-nu, s2-rho, nu-rho
      defval('ocalc','A')
      switch ocalc
        case 'A'
          % What are we staying close to?
          closeto=th0;
          %closeto=mobs;
          fname=fullfile(getenv('IFILES'),'HASHES',...
                         hash([struct2array(orderfields(p)) closeto],'SHA-1'));
          if ~exist(sprintf('%s.mat',fname),'file')
              % Fish for a data vector whose parameter estimates are within X% relative
              % distance of a certain target (this will take a few iterations)
              thhat=NaN; scl=NaN;
              while any(abs(thhat.*scl-closeto)./closeto>0.01) | isnan(thhat)
                  % Simulate from the mean observation
                  p.blurs=Inf; Hx=simulosl(closeto,p);
                  % Calculate the MLE for the data vector - no more gratuitous output
                  [thhat,~,lpars,scl,~,~,Hk]=mleosl(Hx,[],p,[],[],[],[],0);
              end
              [Lgrid,Lcon,thR,xcon,ycon]=...
                  mlelcontosl(Hk,thhat,scl,p,covcont,thcont,runinpar);
              save(fname,'thhat','scl','p','Hk','Lgrid','Lcon','thR','xcon','ycon')
          else
              load(fname)
          end
        case 'B'
          % Scale any data vector simulated from the truth and calculate 
          % the likelihood surface for this observation with mobs, which assumes the
          % process underlying the data is close to MOBS (as it should be). This is
          % quick, but the likelihood surface will likely not be centered on MOBS
          p.blurs=Inf; Hx=simulosl(th0,p); p.blurs=-1;
          shat=nanstd(Hx(:,1));
          Hk=tospec(Hx./shat,p)/(2*pi);
          scl=10.^round(log10(abs(mobs))); scl(1)=shat^2;
          [Lgrid,Lcon,thR,xcon,ycon]=...
              mlelcontosl(Hk,mobs./scl,scl,p,covcont,thcont,runinpar);
        case 'C'
          % Simulate N data vectors from the truth, calculate the 
          % likelihood surface from each observation and their MLE solution, and
          % plot the contours of the average surface 
          N=25;
          for cnd=1:N
              thhat=NaN; 
              while isnan(thhat)
                  p.blurs=Inf; Hx=simulosl(th0,p); p.blurs=-1;
                  [thhat,~,lpars,scl,~,~,Hk]=mleosl(Hx,[],p);
              end
              [Lgridc{cnd},Lconc{cnd},~,xcon,ycon]=...
                  mlelcontosl(Hk,thhat,scl,p,covcont,thcont,runinpar);
          end
          Lgrid=mean(cat(4,Lgridc{:}),4);
          Lcon =mean(cat(4,Lconc{:}), 4);
      end
      
      % Plot the contours; set a grayscale coloring scheme
      figure(2)
      % Remove the points and the cross
      delete([ps o1 o2])
      % Remove the marker at the center indicating the mean observation
      delete(ms)
      % Remove the other annotations
      delete(ttt)
      delete(ep)
      delete(ec)

      % Just reusing scl isn't good enough, now it's consitent with MLEPLOS axes
      sclth0=10.^round(log10(abs(th0)));
      
      % Prepare to plot the loglihood contours on the existing axes
      % Find the pairwise combinations for the cross-plot convention: s2-nu,
      % s2-rho, nu-rho
      pcomb=nchoosek(1:np,2);
      for ind=1:np
          xi=pcomb(ind,1); yi=pcomb(ind,2);
          axes(ah(ind))
          hold on
          % Contour or controuf
          [cont,ch(ind)]=contourf(xcon{ind}/sclth0(xi),ycon{ind}/sclth0(yi),...
                                  Lgrid(:,:,ind),Lcon(ind,:));
          % This doesn't seem to do much now
          ch(ind).LineWidth=0.1;
          % This is slightly better says Olivia
          ch(ind).EdgeAlpha=0.4;
          % set(ch(ind),'EdgeColor',grey(3),'LineStyle','-');
          cmap=flipud(gray(10));
          % Cut off the last bit
          cmap(end,:)=[1 1 1]; 
          colormap(cmap)
          bottom(ch(ind),ah(ind))
          box on
      end
      hold off
      keyboard
      figure(2); figna=figdisp([],sprintf('%s_%s_3','demo2',datum),[],1);
  end
elseif strcmp(Hx,'demo3')
  disp('This does not exist, numbering kept for consistency only')
elseif strcmp(Hx,'demo4')
  defval('thini',[]);
  datum=thini;
  defval('datum',date)
  defval('ifinv',[1 1 1])

  % The number of parameters to solve for
  np=3;

  % Load everything you know about this simulation
  % Looks like more trimming is needed for 'con' rather than 'unc'
  trims=100;
  [th0,thhats,params,covX,~,pix,~,~,obscov,sclcovX,~,covXpix]=osload(datum,trims);

  % Make the plot
  ah=covplos(2,sclcovX,obscov,params,thhats,[],[],'ver',ifinv);

  % Print the figure!
  disp(' ')
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
elseif strcmp(Hx,'demo5')  
  % What th-parameter set? The SECOND argument after the demo id
  defval('thini',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=thini; clear thini
  % What fixed-parameter set? The THIRD argument after the demo id
  defval('params',[]);
  
  % Figure name
  figna=sprintf('%s_%s_%s',mfilename,Hx,date);

  % Simulate data, watch the blurring, verify CHOLCHECK inside
  [Hx,th0,p,k,Hk]=simulosl(th0,params,1);
  
  % Initialize, take defaulted inside MLEOSL for now
  thini=[];

  % Perform the optimization, whatever the quality of the result
  [thhat,covFHh,lpars,scl,thini,p,Hks,k]=mleosl(Hx,thini,p);
  matscl=[scl(:)*scl(:)'];

  % We're going to taper just for what follows
  if isinf(p.blurs)
    Tx=gettaper(p,'cosine',0.1);
    % Replicating the treatment applied to Hk in SIMULOSL
    Hx=Tx(:).*Hx(:,1); Hx=Hx-mean(Hx);
    p.taper=Tx;
    Hk=tospec(Hx,p)/(2*pi);
    if size(Tx)~=1
        Hk=Hk/sqrt(sum(Tx(:).^2))*sqrt(prod(p.NyNx));
    end
    % Note we won't look at the solution, but we use the taper for the residual
  end

  if any(isnan(k(:))); return; end

  % Fisher and Fisher-derived covariance at the truth
  [F0,covF0]=fishiosl(k,th0,params);
  % Fisher and Fisher-derived covariance at the estimate
  % covF=covFHh{1};
   
  % Scaled covariances based on the analytical Hessian at the estimate
  predcov=covFHh{2}./[scl(:)*scl(:)'];
  % Scaled covariances based on the numerical Hessian at the estimate
  obscov=covFHh{3}./[scl(:)*scl(:)'];

  % Quick status report, but note that you get more information in demo4
  disp(sprintf('%s\n',repmat('-',1,97)))
  disp('Analytical and numerical scaled Hessian standard deviations and their ratio')
  disp('   s2      nu     rho')
  disp(sprintf([repmat('%6.3f  ',1,length(obscov)) '\n'],...
                [sqrt(diag(predcov))' ; sqrt(diag(obscov))' ; ...
         	sqrt(diag(predcov))'./sqrt(diag(obscov))']'))
  disp(repmat('-',1,97))

  % Talk again!
  [str0,str2]=osdisp(th0,p);
  % Update like inside MLEOSL
  str2='%13.3g %13.3g %13.5g';
  disp(sprintf(sprintf('%s : %s ',str0,str2),...
	       'Estimated theta',thhat.*scl))
  disp(repmat('-',1,97))

  % Subvert OSDISP for this one example
  osdisp(th0,thhat,1,trilos(lpars{3}),trilos(F0.*matscl),covFHh{3})
  
  % Quick plot, but see also EGGERS3
  clf
  ah=krijetem(subnum(2,3)); delete(ah(4:6)); ah=ah(1:3);

  if ~isinf(p.blurs)
      % Time to rerun LOGLIOS one last time at the solution
      sclh=scl; sclh(1)=1;
      % We had this already, just making sure it checks out if we hadn't changed it
      [L,~,~,momx]=logliosl(k,sclh.*thhat,p,Hks,1);
      diferm(L,lpars{1})
  end

  % Makes an attractive plot that can be used as a judgment for fit
  [cb,xl,t,tt]=mlechiplos(4,Hk,thhat,scl,p,ah,0,th0,covFHh{3});

  disp('FJS here fits also MLECHIPSDOSL')
  disp('FJS here fits the MLELCONTOSL')

  % Print the figure!
  disp(' ')
  figna=figdisp(figna,[],[],1)
elseif strcmp(Hx,'demo6')
    % Simulate something
    [Hx,th0,params]=simulosl;
    % Optimize the smoothness
    tic; [thhat,~,~,scl]=mleosl(Hx,[],params,[],[],[],[1 1 1]); toc
    % Do not optimize the smoothness
    thini=thhat.*scl; thini(2)=th0(2);
    tic; [thhat2,~,~,scl2]=mleosl(Hx,thini,params,[],[],[],[1 0 1]); toc
    disp(sprintf('3-parameter inversion estimate: [%0.2f %0.2f %0.2f]',...
        thhat.*scl))
    disp(sprintf('2-parameter inversion estimate (nu fixed): [%0.2f %0.2f %0.2f]',...
        thhat2.*scl2))
end
