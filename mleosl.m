function varargout=mleosl(Hx,thini,params,algo,bounds,aguess,xver)
% [thhat,covF,logli,thini,scl,params,eflag,oput,grd,hes,Hk,k,options,bounds]=...
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
% algo     'unc' uses FMINUNC for unconstrained optimization
%          'con' uses FMINCON with positivity constraints [default]
%          'klose' simply closes out a run that got stuck
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
% covF     A covariance estimate, from the Fisher matrix, at the estimate
% logli    The maximized value of the likelihood
% thini    The scaled starting guess used in the optimization
% scl      The scaling applied as part of the optimization procedure
% params   The known constants used inside, see above
% eflag    The exit flag of the optimization procedure [bad if 0]
% oput     The output structure of the optimization procedure
% grd      The scaled gradient of the misfit function at the estimate
% hes      The scaled Hessian of the misfit function at the estimate
% Hk       The spectral-domain interface topography after deconvolution
% k        The wavenumbers on which the estimate is actually based
% options  The options used by the optimization procedure
% bounds   The bounds used by the optimization procedure
%
% NOTE: 
%
% A program like EGGERS5 runs 'demo1' in an SPMD loop. Files are opened in
% append mode, except "thzro", which only reflects one lab in that case.
% Writing wires could get cross-checked, messing up the "diagn" files.
%
% EXAMPLE:
%
% p.quart=0; p.blurs=0; p.kiso=NaN;
% [Hx,th0,p]=simulosl([],p,1);
% [thh,covF,L,thi,scl,pp]=mleosl(Hx,[],p,[],[],[],1);
%
% You can stick in partial structures, e.g. only specifying params.kiso
%% Perform a series of N simulations centered on th0
% mleosl('demo1',N,th0,params)
%% Statistical study of a series of simulations done using 'demo1'
% mleosl('demo2','02-Oct-2014') % You can ask for output 
%% Covariance study of a series of simulations
% mleosl('demo4','02-Oct-2014')
%% One simulation and a chi-squared plot
% mleosl('demo5',th0,params)
%
% Last modified by fjsimons-at-alum.mit.edu, 10/14/2016

% NEED TO CHANGE THE k(~~k) to proper accounting for kiso

if ~isstr(Hx)
  defval('algo','unc')
  % The necessary strings for formatting FJS see OSDISP
  str0='%27s';
  str1='%12.0e ';
  str2='%12.5g ';

  % Supply the needed parameters, keep the givens, extract to variables
  fields={               'dydx','NyNx','blurs','kiso','quart'};
  defstruct('params',fields,...
	    {                      [20 20]*1e3,sqrt(length(Hx))*[1 1],2,NaN,0});
  struct2var(params)

  % These bounds are physically motivated...
  if strcmp(algo,'con')
    % Parameters for FMINCON in case that's what's being used, which is recommended
    defval('bounds',{[],[],... % Linear inequalities
                     [],[],... % Linear equalities
                     [var(Hx)/100 0.15  sqrt(prod(dydx))],... % Lower bounds
                     [var(Hx)*100 8.00  max(2.5e5,min(dydx.*NyNx))],... % Upper bounds
                     []}); % Nonlinear (in)equalities
    % nu bound: 
    disp(sprintf('NU lower bound %5.2f',bounds{5}(2)))
    disp(sprintf('NU upper bound %5.2f',bounds{6}(2)))
  else
    bounds=[];
  end
  
  % Being extra careful or not?
  defval('xver',0)

  % The parameters used in the simulation for demos, or upon which to base "thini"
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
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,size(scl))),...
		 'Scaling',scl))
  end
  disp(sprintf(sprintf('%s : %s',str0,repmat(str2,size(thini))),...
	       'Starting theta',thini))

  % Now scale so the minimization doesn't get into trouble
  thini=thini./scl;
  
  defval('taper',0)
  if taper==1
    % Were going to want to make a 2D taper - any taper
    disp(sprintf('%s with TAPERING, DO NOT DO THIS YET',upper(mfilename)))
    NW=2;
    E1=dpss(NyNx(1),NW,1);
    E2=dpss(NyNx(2),NW,1);
    Tx=repmat(E1,1,NyNx(2)).*repmat(E2',NyNx(1),1);
    % But should still watch out for the spectral gain I suppose, this isn't
    % done just yet, although it appears already properly normalized
    % However, this looks better now, doesn't it?
    Tx=Tx*sqrt(prod(NyNx));
    % Not doing anything still amounts to saying Tx=1
  else
    Tx=1;
  end

  % Create the appropriate wavenumber axis
  k=knums(params);
  
  % We get into the habit of never involving the zero-wavenumber
  knz=(~~k);

  % Always demean the data sets
  Hx(:,1)=Hx(:,1)-mean(Hx(:,1));

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

  if xver==1 && blurs>-1 && blurs<2
    % Using the analytical gradient in the optimization is not generally a good
    % idea but if the likelihoods aren't blurred, you can set this option to
    % 'on' and then you can verify that the numerical calculations match the
    % analytics. According to the manual, "solvers check the match at a
    % point that is a small random perturbation of the initial point". My
    % own disp output (further below) provides comparisons at the estimate.
    options.GradObj='off';
    % When it's "on", you can perform the derivative check
    if strcmp(options.GradObj,'on');
      % The "if" statement wasn't strictly necessary, since it wouldn't actually
      % check the derivatives is GradObj wasn't on, even when LOGLIOS has them.
      options.DerivativeCheck='on';
    end
  end

  % And find the MLE! Work on scaled parameters
  try
    switch algo
     case 'unc'
      % disp('Using FMINUNC for unconstrained optimization of LOGLIOSL')
      t0=clock;
      [thhat,logli,eflag,oput,grd,hes]=...
	  fminunc(@(theta) logliosl(theta,params,Hk,k,scl),...
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
            fmincon(@(theta) logliosl(theta,params,Hk,k,scl),...
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
       varargout=cellnan(nargout,1,1);
       varargout{end-1}=options;
       varargout{end}=bounds;
       return
    end
  catch
    % If something went wrong, exit gracefully
    varargout=cellnan(nargout,1,1);
    return
  end

  % Parameter covariance as calculated from numerical Hessian at estimate
  covh=hes2cov(-hes,scl,length(k(knz))/2);

  % Parameter covariance as calculated from unblurred Fisher matrix at estimate
  % I suppose, in the single-variable case, we could produce blurred versions
  [covF,F]=covthosl(thhat,k(knz),scl); 
    
  % Analytical calculations of the Hessian neglect blurring, and thus,
  % are expected to be close to numerical results only without blurring
  % I don't yet understand the behavior if we DO compare them - the
  % Hessians are off, but the analytic Fishers remain right on; maybe the
  % blurring removes the variability such that numerical blurring ends up
  % producing results that are closer to the analytical Fishers
  if xver==1 && blurs>-1 && blurs<2
    % The number of parameters that are being solved for
    np=length(thhat);
    % The number of unique entries in an np*np symmetric matrix
    npp=np*(np+1)/2;

    % Analytic (unblurred) Hessian, scaled for numerical comparison
    H=Hessiosl(k(knz),thhat.*scl,Hk(knz)).*[scl(:)*scl(:)'];

    % Parameter covariance as calculated from analytical Hessian at estimate
    covH=hes2cov(H,scl,length(k(knz))/2);
    
    % Analytic (unblurred) gradient, scaled for numerical comparison
    gros=-nanmean(gammakosl(k(knz),thhat.*scl,params,Hk(knz)))'.*scl(:);

    % Compare the analytic Hessian with the numerical Hessian and with
    % the Hessian expectation, which is the Fisher, at the estimate, and
    % compare the analytic gradient with the numerical gradient
    if all(thhat>0)
      str3=repmat('%16g ',1,npp);
      str4=repmat('%16g ',1,np);
      disp(sprintf('\n%16s\n','At the ESTIMATE:'));
      disp(sprintf(sprintf('Log likelihood    %s',str3),logli))
      disp(' ')
      disp(sprintf(sprintf('Numericl Gradient %s',str4),grd))
      disp(sprintf(sprintf('Analytic Gradient %s',str4),gros))
      disp(' ')
      disp(sprintf(sprintf('Numerical Hessian %s',str3),trilos(hes)))
      disp(sprintf(sprintf('Analyticl Hessian %s',str3),trilos(-H )))
      disp(sprintf(sprintf('Analytical Fisher %s',str3),trilos(F  )))
      disp(' ')
      disp(sprintf(sprintf('Cov (Numer Hess.) %s',str3),trilos(covh)))
      disp(sprintf(sprintf('Cov (Analy Hess.) %s',str3),trilos(covH)))
      disp(sprintf(sprintf('Cov (Analy Fish.) %s',str3),trilos(covF)))
    end
  end

  % Talk!
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Numer Hessi std',sqrt(diag(covh))))
  disp(sprintf(sprintf('%s : %s\n ',str0,repmat(str2,size(thhat))),...
	       'Anal Fisher std',sqrt(diag(covF))))
  if xver==1
    disp(sprintf('%8.3gs per %i iterations or %8.3gs per %i function counts',...
                 ts/oput.iterations*100,100,ts/oput.funcCount*1000,1000))
  else
    disp(sprintf('\n'))
  end

  % Generate output as needed
  varns={thhat,covF,logli,thini,scl,params,eflag,oput,grd,hes,Hk,k,options,bounds};
  varargout=varns(1:nargout);
elseif strcmp(Hx,'demo1')
  % Runs a series of simulations. See 'demo2' to display them.
  % If you run this again on the same date, the files THINI and
  % THHAT get appended, but a blank THZERO is created. 
  defval('thini',[]);
  % How many simulations? The SECOND argument, after the demo id.
  N=thini; clear thini
  defval('N',500)
  more off
  % What th-parameter set? The THIRD argument, after the demo id
  defval('params',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=params; clear params
  % What fixed-parameter set? The FOURTH argument, after the demo id
  defval('algo',[]);
  % If there is no preference, then that's OK, it gets taken care of
  params=algo; clear algo
  % What algorithm? The FIFTH argument, after the demo id
  defval('bounds',[]);
  % If there is no preference, then that's OK, it gets taken care of
  algo=bounds; clear bounds
    
  % The number of parameters to solve for
  np=3;
  
  % Open files and return format strings
  [fid0,fid1,fid2,fid3,fmt1,fmt2,fmt3,fmtf,fmte,fmtd,fmtc]=...
      osopen(np);

  % Do it!
  good=0; 
  % Initialize the average Hessian
  avH=zeros(np,np);

  % Set N to zero to simply close THZERO out
  for index=1:N
    % Simulate data from the same lithosphere, watch the blurring
    [Hx,th0,p,k,Hk]=simulosl(th0,params);

    % Check the dimensions of space and spectrum are right
    difer(length(Hx)-length(k(:)),[],[],NaN)

    % Form the maximum-likelihood estimate, pass on the params, use th0
    % as the basis for the perturbed initial values. Remember hes is scaled.
    t0=clock;
    [thhat,covF,logli,thini,scl,p,e,o,gr,hs,Hk,k,ops,bnds]=...
	mleosl(Hx,[],p,algo,[],th0);
    ts=etime(clock,t0);

    % Initialize the THZRO file... note that the bounds may change
    % between simulations, and only one gets recorded here
    if index==1 && labindex==1
      oswzerob(fid0,th0,p,ops,bnds,fmt1,fmt2)
    end
    
    % If a model was found, keep the results, if not, they're all NaNs
    % Ignore the fact that it may be at the maximum number of iterations
    % e=1

    % IF NUMBER OF FUNCTION ITERATIONS IS TOO LOW DEFINITELY BAD
    itmin=0;
    % A measure of first-order optimality (which in this
    % unconstrained case is the infinity norm of the gradient at the
    % solution)  
    % FJS to update what it means to be good - should be in function of
    % the data size as more precision will be needed to navigate things
    % with smaller variance! At any rate, you want this not too low.
    optmin=Inf;
    % Maybe just print it and decide later? No longer e>0 as a condition.
    % e in many times is 0 even though the solution was clearly found, in
    % other words, this condition IS a way of stopping with the solution
    % Remember that the correlation coefficient can be negative or zero!
    % The HS is not always real, might be all the way from the BLUROS?
    % Because if there are NaNs or not estimate it won't work
    try
      % Maybe I'm too restrictive in throwing these out? Maybe the
      % Hessian can be slightly imaginary and I could still find thhat
      if isreal([logli gr(:)']) ...
	    && all(thhat>0) ...
	    && all(~isnan(thhat)) ...
	    && o.iterations > itmin ...
	    && o.firstorderopt < optmin
	good=good+1;
	% Build the average of the Hessians for printout later
	avH=avH+hs./[scl(:)*scl(:)'];
	% Reapply the scaling before writing it out
	fprintf(fid1,fmt1,thhat.*scl);
	fprintf(fid2,fmt1,thini.*scl);
	% Here we compute and write out the moments of the Xk
	[L,~,momx]=logliosl(thhat,p,Hk,k,scl);
	% Print the optimization results and diagnostics to a different file 
	oswdiag(fid3,fmt1,fmt3,logli,gr,hs,thhat,thini,scl,ts,e,o,....
		var(Hx),momx,covF)
      end
    end
  end
  % If there was any success at all, finalize the THZERO file 
  % If for some reason this didn't end well, do an N==0 run.

  % Initialize if all you want is to close the file
  if N==0
    [Hx,th0,p,k]=simulosl(th0,params); 
    good=1; avH=avH+1; 
    [~,~,~,~,~,~,~,~,~,~,~,~,ops,bnds]=mleosl(Hx,[],[],'klose');
    oswzerob(fid0,th0,p,ops,bnds,fmt1,fmt2)
  end
  
  if good>=1 
    % This is the scaling based on the truth which we use here 
    sclth0=10.^round(log10(th0));

    % This is the average of the Hessians, should be closer to the Fisher
    avH=avH.*[sclth0(:)*sclth0(:)']/good;

    % Now compute the theoretical covariance and scaled Fisher, avoid zero-k
    [covF,F]=covthosl(th0./sclth0,k(knz),sclth0);
    % Of course when we don't have the truth we'll build the covariance
    % from the single estimate that we have just obtained. This
    % covariance would then be the only thing we'd have to save.
    if labindex==1
      oswzeroe(fid0,sclth0,avH,good,F,covF,fmtc,fmte,fmtf)
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

  % Load everything you know about this simulation
  [th0,thhats,p,covF,covHav,covHts,~,~,~,~,momx,covthpix]=osload(datum);

  % Report the findings of all of the moment parameters
  disp(sprintf('m(m(Xk)) %f m(v(Xk)) %f\nm(magic) %f v(magic) %f',...
	      mean(momx),var(momx(:,end))))
  
  % Plot it all - perhaps some outlier selection?
  disp(sprintf('\n'))
  figure(1)
  fig2print(gcf,'landscape')
  clf
  trims=98;
  % disp(sprintf('%s estimates trimmed at %i percentile',...
  %      upper(mfilename),trims))
  
  % covF is the Fisher-based covariance evaluate at the truth
  % covthpix is the Fisher-based covariance at one of the randomly picked estimates

  % If the above two are close, we need to start using the second one
  % Yep!
  [ah,ha]=mleplos(trimit(thhats,trims,1),th0,covF,covHav,covthpix,[],[],p,...
                  sprintf('MLEOSL-%s',datum));
  
  % Soon, we will be able to start loading different numbers

  % Return some output, how about just the empirically observed
  % means and covariance matrix of the estimates, and the number of
  % reported trials
  nobs=size(thhats,1);
  mobs=mean(thhats,1);
  cobs=cov(thhats);
  varns={cobs,mobs,nobs,th0,p,momx};
  varargout=varns(1:nargout);

  % Print the figure! Don't forget the degs.pl script
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
  
  % Being extra careful or not?
  defval('xver',0)

  if xver==1
    % Take a look a the distribution of the residual moments
    % See RB X, p. 51 about the skewness of a chi-squared - just sayin'.
    % We don't change the number of degrees of freedom! If you have used
    % twice the number, and given half correlated variables, you do not
    % change the variance, that is the whole point. Unlike in COVTHOSL
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
  [th0,thhats,params,covF,~,~,E,v,obscov,sclcov]=osload(datum);

  % Make the plot
  ah=covplos(2,sclcov,obscov,covF,params,thhats,th0,[],[],'ver');

  % Make the plot
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
elseif strcmp(Hx,'demo5')  
  % What th-parameter set? The SECOND argument after the demo id
  defval('thini',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=thini; clear thini
  % What fixed-parameter set? The THIRD argument after the demo id
  defval('params',[]);
  
  % Figure name
  figna=sprintf('%s_%s_%s',mfilename,Hx,date);
  
  % Simulate data, watch the blurring, verify COLCHECK inside
  [Hx,th0,p,k,Hk]=simulosl(th0,params,1);
  
  % Initialize, take defaulted inside MLEOSL for now
  thini=[];

  % Perform the optimization, whatever the quality of the result
  [thhat,~,logli,thini,scl,p,e,o,gr,hs]=mleosl(Hx,thini,p);

  % Take a look at the unblurred gradient purely for fun, they should be
  % so small as to be immaterial
  grobs=-nanmean(gammakosl(k,thhat.*scl,p,Hk))';
  
  % Take a look at the unblurred theoretical covariance at the estimate,
  % to compare to the observed blurred Hessian; in the other demos we
  % compare how well this works after averaging
  [covthat,F]=covthosl(thhat,k,scl);

  % Collect the theoretical covariance for the truth for the title
  covth=covthosl(th0./scl,k,scl);
 
  % Take a look at the scaled Fisher to compare with the scaled Hessian  
  F;
  hs;
  grobs;
  
  % Take a look at the scaled covariances
  predcov=covthat./[scl(:)*scl(:)'];
  % And compare to the inverse of the scaled Hessians in the full plane
  % Watch as this will change with the anticipated halfplane changes
  obscov=inv(hs)/length(k(:))*2;
  % These two should compare reasonably well in the unblurred case, I
  % would have thought - but of course it's ignoring stochastic
  % variability. If we can use the Hessian for the uncertainty estimation
  % we can do this for the cases where we can't come up with a
  % theoretical covariance, not even an unblurred one. Check std's.
  % Maybe should do this a hundred times?
  disp(sprintf('%s\n',repmat('-',1,97)))
  disp('predicted and observed scaled standard deviations and their ratio')
  disp(sprintf([repmat('%6.3f  ',1,length(obscov)) '\n'],...
	       [sqrt(diag(predcov))' ; sqrt(diag(obscov))' ; ...
		sqrt(diag(predcov))'./sqrt(diag(obscov))']'))
  disp(repmat('-',1,97))
  % Talk again!
  [str0,str2]=osdisp(th0,p);
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(repmat('-',1,97))
  
  % Quick plot, but see OLHEDESIMONS5
  clf
  ah=krijetem(subnum(2,3)); delete(ah(4:6)); ah=ah(1:3);

  % Maybe we should show different covariances than the predicted ones??

  % Time to rerun LOGLIOS one last time at the solution
  [L,~,momx]=logliosl(thhat,p,Hk,k,scl);

  disp('fjs feed this into the next one')
  disp('fjs here fits the likelihood contours')

  % Better feed this to the next code, now it's redone inside
  mlechiplos(4,Hk,thhat,scl,p,ah,0,th0,covth);
   
  % Better put in mlelcontsol and mlepsdosl

  % Print the figure! Don't forget the degs.pl script
  figna=figdisp(figna,[],[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
end 
