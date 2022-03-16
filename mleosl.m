function varargout=mleosl(Hx,thini,params,algo,bounds,aguess,xver)
% [thhat,covFHh,lpars,scl,thini,params,Hk,k]=...
%          MLEOSL(Hx,thini,params,algo,bounds,aguess,xver)
%
% Performs a maximum-likelihood estimation for SINGLE FIELDS as in
% Olhede & Simons (2013) by minimization using FMINUNC/FMINCON.
%
% INPUT:
%
% Hx       Real-valued column vector of unwrapped spatial-domain quantities 
% thini    An unscaled starting guess for the parameter vector with elements:
%          [          s2 nu rho], see SIMULOSL. If you leave this value
%          blank, then you will work from the perturbed "aguess" 
% params   A parameter structure with constants assumed known, see SIMULOSL
%          [dydx NyNx blurs kiso] in the units of 
%          m (2x), "nothing" (3x), rad/m, "nothing", namely, in order:
%          blurs  0 Don't blur likelihood using the Fejer window
%                 N Blur likelihood using the [default: N=2] resampled Fejer window
%                -1 Blur likelihood using the exact BLUROSY procedure
%          kiso   wavenumber beyond which we are not considering the likelihood
%          quart 1 quadruple, then QUARTER the spatial size
%                0 size as is, watch for periodic correlation behavior
%          taper 1 a certain taper
%                0 no taper
% algo     'unc' uses FMINUNC for unconstrained optimization
%          'con' uses FMINCON with positivity constraints [default]
%          'klose' simply closes out a run that got stuck [defaulted when needed ]
% bounds    A cell array with those positivity constraints [defaulted]
% aguess    A parameter vector [s2 nu rho] that will be used in
%           simulations for demo purposes, and on which "thini" will be
%           based if that was left blank. If "aguess" is blank, there is
%           a default. If "thini" is set, there is no need for "aguess"
% xver      Conduct extra verification steps
%
% OUTPUT:
%
% thhat    The maximum-likelihood estimate of the vector [scaled]:
%          [s2 nu rho], in units of variance, "nothing", and distance, see SIMULOSL
% covFHh   The covariance estimates:
%          covFHh{1} from Fisher matrix AT the estimate [FISHIOSL]
%          covFHh{2} from analytical Hessian matrix AT the estimate [HESSIOSL]
%          covFHh{3} from numerical Hessian matrix NEAR the estimate [FMINUNC/FMINCON]
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
% p.quart=0; p.blurs=Inf; p.kiso=NaN; clc; [Hx,th,p]=simulosl([],p,1); mleosl(Hx,[],p,[],[],[],1);
%
% You can stick in partial structures, e.g. only specifying params.kiso
%
%% Perform a series of N simulations centered on th0
% mleosl('demo1',N,th0,params)
%
%% Statistical study of a series of simulations using MLEPLOS
% mleosl('demo2','20-Jun-2018')
%
%% Covariance study of a series of simulations using COVPLOS
% mleosl('demo4','20-Jun-2018')
%
%% One simulation and a chi-squared plot using MLECHIPLOS
% mleosl('demo5',th0,params)
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/10/2022

if ~isstr(Hx)
  defval('algo','unc')
  % The necessary strings for formatting, see OSDISP and OSANSW
  str0='%18s';
  str1='%13.0e ';
  str2='%13.5g %13.5g %13.5g';
  str2='%13.0f %13.2f %13.0f';
  str3s='%13s ';

  % Supply the needed parameters, keep the givens, extract to variables
  fields={               'dydx','NyNx','blurs','kiso','quart','taper'};
  defstruct('params',fields,...
	    {                      [20 20]*1e3,sqrt(length(Hx))*[1 1],Inf,NaN,0,1});
  struct2var(params)

  % These bounds are physically motivated...
  if strcmp(algo,'con')
    % Parameters for FMINCON in case that's what's being used, which is recommended
    defval('bounds',{[],[],... % Linear inequalities
                     [],[],... % Linear equalities
                     [var(Hx)/100 0.15  sqrt(prod(dydx))],... % Lower bounds
                     [var(Hx)*100 8.00  max(2.5e5,min(dydx.*NyNx))],... % Upper bounds
                     []}); % Nonlinear (in)equalities
  else
    bounds=[];
  end
  
  % Being extra careful or not?
  defval('xver',1)

  % The parameters used in the simulation for demos, or upon which to base "thini"
  % Check Vanmarcke 1st edition for suggestions on initial rho
  defval('aguess',[var(Hx) 2.0 sqrt(prod(dydx.*NyNx))/5]);
  % Scale the parameters by this factor; fix it unless "thini" is supplied
  defval('scl',10.^round(log10(abs(aguess))));

  % Unless you supply an initial value, construct one from "aguess" by perturbation
  nperturb=0.25;
  % So not all the initialization points are the same!!
  
  defval('thini',abs((1+nperturb*randn(size(aguess))).*aguess))
  % If you brought in your own initial guess, need an appropriate new scale
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
  
  % Now scale so the minimization doesn't get into trouble
  thini=thini./scl;
  
  % I do believe that when params.blurs=Inf we should be tapering since
  % then the data are no longer periodic on the space grid...
  if taper==1
    % We're going to want to make a 2D taper - any taper
    % disp(sprintf('%s with TAPERING',upper(mfilename)))
    Tx=ones(params.NyNx);
    %NW=2;
    %E1=dpss(NyNx(1),NW,1); E2=dpss(NyNx(2),NW,1);
    %E1=shanning(NyNx(1),0.2); E2=shanning(NyNx(2),0.2); E1=E1/sqrt(E1'*E1); E2=E2/sqrt(E2'*E2);
    %Tx=repmat(E1,1,NyNx(2)).*repmat(E2',NyNx(1),1);
    % But should still watch out for the spectral gain I suppose, this isn't
    % done just yet, although it appears already properly normalized
    % However, this looks better now, doesn't it?
    %Tx=Tx*sqrt(prod(NyNx));
    % Just the rim removal
    Tx(1,:)=0; Tx(end,:)=0; Tx(:,1)=0; Tx(:,end)=0;
    
    % A whole Tukey/Hanning window
    % [w,wl,wr]=fhanning(floor(sqrt(prod(params.NyNx)))/5);
    % Tx(:,1:length(wl))=repmat(wl',size(Tx,1),1);
    % Tx(:,end-length(wl)+1:end)=repmat(wr',size(Tx,1),1);
    % Tx(1:length(wr),:)=Tx(1:length(wr),:).*repmat(wl,1,size(Tx,2));
    % Tx(end-length(wr)+1:end,:)=Tx(end-length(wr)+1:end,:).*repmat(wr,1,size(Tx,2));
    
    % Make sure to renormalize the taper
    Tx=Tx*sqrt(prod(NyNx)/[Tx(:)'*Tx(:)]);
  else
    Tx=1;
  end

  % Create the appropriate wavenumber axis
  k=knums(params);
  
  % We could get into the habit of never involving the zero-wavenumber
  knz=(~~k);

  % Always demean the data sets
  Hx(:,1)=Hx(:,1)-mean(Hx(:,1));
  % FJS think about deplaning as well

  % Turn the observation vector to the spectral domain
  % Watch the 2pi in SIMULOSL
  Hk(:,1)=tospec(Tx(:).*Hx(:,1),params)/(2*pi);

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
      t0=clock;
      [thhat,logli,eflag,oput,grd,hes]=...
	  fminunc(@(theta) logliosl(k,theta,scl,params,Hk,xver),...
		  thini,options);
      ts=etime(clock,t0);
     case 'con'
      % New for FMINCON
      options.Algorithm='active-set';
      % disp('Using FMINCON for constrained optimization of LOGLIOSL')
      t0=clock;
      % See M. K. Stein p. 173 when differentiability parameter maxes out
      % Also check for when this crucial parameter hits the lower bound
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
        [thhat,logli,eflag,oput,lmd,grd,hes]=...
            fmincon(@(theta) logliosl(k,theta,scl,params,Hk,xver),...
                    thini,...
                    bounds{1},bounds{2},bounds{3},bounds{4},...
                    bounds{5}./scl,bounds{6}./scl,bounds{7},...
                    options);
        % Try resetting the offending parameter nu by a serious kick
        thini(2)=thini(2)/[1+1/4-rand/2];
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
    varargout=cellnan(nargout,1,1);
    return
  end

  % It is not impossible that a solution is reached which yields a
  % negative rho - which only appears in the square in MATERNOS. But if
  % we're going to calculate (approximately blurred) analytical
  % gradients and Hessians (even using exact blurring of the spectral
  % densities) we are going to be using MATERNOSY, which will complain...
  if thhat(1)<0
    error(sprintf('%s Negative variance',upper(mfilename)))
  end
  if thhat(2)<0
    thhat(2)=-thhat(2);
  end
  if thhat(3)<0
    thhat(3)=-thhat(3);
  end

  % Degrees of freedom for full-wavenumber domain (redundant for real data)
  % Not including the zero wavenumber, since LOGLIOS doesn't either
  df=length(k(~~k))/2; 
  matscl=[scl(:)*scl(:)'];

  % Watch out for singularity or scaling warnings, they are prone to pop up
  
  % Covariance from FMINUNC/FMINCON's numerical scaled Hessian NEAR estimate
  covh=inv(hes./matscl)/df;
  
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
  [F,covF]=fishiosl(k,thhat.*scl,xver);

  % Analytic (poorly blurred) Hessian AT the estimate, and derived covariance
  [H,covH]=hessiosl(k,thhat.*scl,params,Hk,xver);

  % FJS how about a step further, use F-1 H F-T to get any influence at all
  covFHF=inv(F)*[-H]*inv(F)/df;
  covFhF=inv(F)*[hes./matscl]*inv(F)/df;

  % Analytical calculations of the gradient and the Hessian poorly represent
  % the blurring (though it's much better than not trying at all), and thus,
  % are expected to be close to numerical results only without blurring
  if xver==1 
    % Analytic (poorly blurred) gradient, scaled for numerical comparison
    gros=gammiosl(k,thhat.*scl,params,Hk,xver).*scl(:);

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
    disp(sprintf(sprintf(' Numerical Hessian : %s',str3),trilos(hes)))
    disp(sprintf(sprintf(' Analyticl Hessian : %s',str3),trilos(-H.*matscl)))
    disp(sprintf(sprintf(' Analytical Fisher : %s',str3),trilos( F.*matscl)))

    disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3s,1,npp)),...
	       ' ','C(s2,s2)','C(nu,nu)','C(rho,rho)','C(s2,nu)','C(s2,rho)','C(nu,rho)'))
    disp(sprintf(sprintf(' Cov (Numer Hess.) : %s',str3),trilos(covh)))
    disp(sprintf(sprintf(' Cov (Analy Hess.) : %s',str3),trilos(covH)))
    disp(sprintf(sprintf(' Cov (Analy Fish.) : %s',str3),trilos(covF)))
    disp(sprintf(sprintf(' Cov ( FishHFish.) : %s',str3),trilos(covFHF)))
    disp(sprintf(sprintf(' Cov ( FishHFish.) : %s',str3),trilos(covFhF)))
    disp(sprintf('%s',repmat('_',119,1)))
  end

  % Talk!
  disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3s,1,np)),...
	       ' ','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str2),...
	       'Estimated theta',thhat.*scl))
  disp(' ')
  disp(sprintf(sprintf('%s : %s ',str0,str2),...
	       'Numer Hessi std',sqrt(diag(covh))))
  disp(sprintf(sprintf('%s : %s ',str0,str2),...
	       'Analy Hessi std',sqrt(diag(covH))))
  disp(sprintf(sprintf('%s : %s\n ',str0,str2),...
	       'Anal Fisher std',sqrt(diag(covF))))
  disp(sprintf(sprintf('%s : %s ',str0,str2),...
	       ' FishHFish. std',sqrt(diag(covFHF))))
  disp(sprintf(sprintf('%s : %s\n ',str0,str2),...
	       ' FishHFish. std',sqrt(diag(covFhF))))
  if xver==1 | xver==0
    disp(sprintf('%s\n',repmat('_',119,1)))
    disp(sprintf('%8.3gs per %i iterations or %8.3gs per %i function counts',...
                 ts/oput.iterations*100,100,ts/oput.funcCount*1000,1000))
  end

  % Here we compute the moment parameters and recheck the likelihood
  [L,~,Hagain,momx,vr]=logliosl(k,thhat,scl,params,Hk,xver);
  diferm(L,logli)
  diferm(Hagain,H)

  % Reorganize the output into cell arrays
  covFHh{1}=covF;
  covFHh{2}=covH;
  covFHh{3}=covh;
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
  varns={thhat,covFHh,lpars,scl,thini,params,Hk,k};
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
  % The SIXTH argument, after the demo id
  defval('aguess',[])
  % If there is no preference, then that's OK, it gets taken care of
  bounds=aguess; clear aguess
  % The SEVENTH argument, after the demo id
  defval('xver',[])
  % If there is no preference, then that's OK, it gets taken care of
  aguess=xver; clear xver

  % You can't stick in an EIGHT argument so you'll have to default 
  defval('xver',0)

  % What you make of all of that if there hasn't been a number specified
  defval('N',500)

  % The number of parameters to solve for
  np=3;
  
  % Open 'thzro', 'thini', 'thhat' and 'diagn' files and return format strings
  [fids,fmts,fmti]=osopen(np);

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

    % Form the maximum-likelihood estimate, pass on the params, use th0
    % as the basis for the perturbed initial values. Remember hes is scaled.
    t0=clock;
    % Put a try catch here to catch the rare negative variance??
    [thhat,covFHh,lpars,scl,thini,p,Hk,k]=mleosl(Hx,[],p,algo,[],th0,xver);
    ts=etime(clock,t0);

    % Initialize the THZRO file... note that the bounds may change
    % between simulations, and only one gets recorded here
    if index==1 && labindex==1
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
	% Build the AVERAGE of the Hessians for printout later
	avhsz=avhsz+lpars{3}./[scl(:)*scl(:)'];
	% Reapply the scaling before writing it out
	fprintf(fids(2),fmts{1},thhat.*scl);
	fprintf(fids(3),fmts{1},thini.*scl);
        % We don't compare the second and third outputs of LOGLIOSL since these are
        % analytical, poorly approximately blurred, derivatives, and we be
        % writing the numerical versions. Be aware that covFHh{3} is the
        % current favorite covariance estimate on the parameters!
	% Print optimization results and diagnostics to different file 
	oswdiag(fids(4),fmts,lpars,thhat,thini,scl,ts,var(Hx),covFHh{3})
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
    % This is the scaling based on the truth which we use here 
    sclth0=10.^round(log10(th0));

    % This is the AVERAGE of the numerical Hessians, should be closer to the Fisher
    avhsz=avhsz.*[sclth0(:)*sclth0(:)']/good;

    % Now compute the Fisher and Fisher-derived covariance at the truth
    [F0,covF0]=fishiosl(k,th0);
    matscl=[sclth0(:)*sclth0(:)'];

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

  % The number of parameters to solve for
  np=3;

  % Looks like more trimming is needed for 'con' rather than 'unc'
  trims=100;
  % Load everything you know about this simulation
  [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(datum,trims);

  % Report the findings of all of the moment parameters
  disp(sprintf('\nm(m(Xk)) %f m(v(Xk)) %f\nm(magic) %f v(magic) %f',...
	      mean(momx),var(momx(:,end))))
  
  % Plot it all - perhaps some outlier selection?
  figure(1)
  fig2print(gcf,'landscape')
  clf

  % If the above two are close, we need to start using the second one
  % Let us feed it the Fisher covariance at the TRUTH
  [ah,ha]=mleplos(thhats,th0,covF0,covavhs,covXpix,[],[],p, ...
                  sprintf('MLEOSL-%s',datum),thpix);
  
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
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],2);
  
  % Being extra careful or not?
  defval('xver',0)

  if xver==1
    % Take a look a the distribution of the residual moments
    % This now is a full part of MLECHIPLOS
    % See RB X, p. 51 about the skewness of a chi-squared - just sayin'.
    % We don't change the number of degrees of freedom! If you have used
    % twice the number, and given half correlated variables, you do not
    % change the variance, that is the whole point. Unlike in FISHIOSL
    % where you make an analytical prediction that does depend on the
    % number and which therefore you need to adjust.
    k=knums(p); varpred=8/[length(k(~~k))];

    figure(2)
    clf
    fig2print(gcf','portrait')
    ahh(1)=subplot(121);
    histfit(momx(:,3)); 
    [m,s]=normfit(momx(:,3));
    disp(sprintf('mean %f predicted mean 1 \nstdv %s predicted stdv %s',m,s,sqrt(varpred)))
    shrink(ahh(1),1,1.5)
    xl(1)=xlabel('histogram of the residual moments');
    ahh(2)=subplot(122);
    qqplot((momx(:,3)-1)/sqrt(varpred)); axis image; grid on; box on
    refline(1,0)
    movev(ahh,-0.1)
    t=ostitle(ahh,p,sprintf('MLEOSL-%s',datum)); movev(t,1)
    % Could also do, as these quantities should be very close of course
    % qqplot(momx(:,2),momx(:,3)); axis image; refline(1,0); grid on
    % Then use NORMTEST to ascertain the veracity... don't bother with the
    % Nyquist wavenumbers, there will be very few, but take out the zero

    % Predicted expected value is one.
    [a,b,c,d]=normtest(momx(:,3),1,varpred);
  end
elseif strcmp(Hx,'demo3')
  disp('This does not exist, numbering kept for consistency only')
elseif strcmp(Hx,'demo4')
  defval('thini',[]);
  datum=thini;
  defval('datum',date)

  % The number of parameters to solve for
  np=3;

  % Load everything you know about this simulation
  % Looks like more trimming is needed for 'con' rather than 'unc'
  trims=100;
  [th0,thhats,params,covX,~,pix,~,~,obscov,sclcovX,~,covXpix]=osload(datum,trims);

  % Make the plot
  ah=covplos(2,sclcovX,obscov,params,thhats,[],[],'ver');

  % Print the figure!
  disp(' ')
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],2);
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
  [thhat,covFHh,lpars,scl,thini,p,Hk,k]=mleosl(Hx,thini,p);
  matscl=[scl(:)*scl(:)'];

  % Fisher and Fisher-derived covariance at the truth
  [F0,covF0]=fishiosl(k,th0);
  % Fisher and Fisher-derived covariance at the estimate
  % covF=covFHh{1};
  % Those two are close of course, and of not much intrinsic interest anymore
   
  % Make sure there isn't a factor of two in-between covFHh{1} and covFHh{2}
  % Sometimes the blurring makes it look like that
  % Should be testing that these are all closer together the larger the
  % data set is 
  
  % Scaled covariances based on the analytical Hessian at the estimate
  predcov=covFHh{2}./[scl(:)*scl(:)'];
  % Scaled covariances based on the numerical Hessian at the estimate
  obscov=covFHh{3}./[scl(:)*scl(:)'];

  % Quick status report, but note that you get more information in demo4
  disp(sprintf('%s\n',repmat('-',1,97)))
  disp('analytical and numerical scaled Hessian standard deviations and their ratio')
  disp(sprintf([repmat('%6.3f  ',1,length(obscov)) '\n'],...
	       [sqrt(diag(predcov))' ; sqrt(diag(obscov))' ; ...
		sqrt(diag(predcov))'./sqrt(diag(obscov))']'))
  disp(repmat('-',1,97))
  % Talk again!
  [str0,str2]=osdisp(th0,p);
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(repmat('-',1,97))

  % Subvert OSDISP for this one example
  osdisp(th0,thhat,1,trilos(lpars{3}),trilos(F0.*matscl),covFHh{3})
  
  % Quick plot, but see also EGGERS3
  clf
  ah=krijetem(subnum(2,3)); delete(ah(4:6)); ah=ah(1:3);

  % Maybe we should show different covariances than the predicted ones??

  % Time to rerun LOGLIOS one last time at the solution
  % Do not collect the analytical gradient and Hessian, since these are
  % not the observed blurred gradients
  [L,~,~,momx]=logliosl(k,thhat,scl,p,Hk,1);
  % We had this already, just making sure it checks out
  diferm(L,lpars{1})

  % Makes an attractive plot that can be used as a judgment for fit
  mlechiplos(4,Hk,thhat,scl,p,ah,0,th0,covFHh{3});
  
  disp('FJS here fits also MLECHIPSDOSL')
  disp('FJS here fits the likelihood contours')

  % Print the figure!
  disp(' ')
  figna=figdisp(figna,[],[],2);
end
