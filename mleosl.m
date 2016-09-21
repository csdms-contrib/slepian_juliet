function varargout=mleosl(Hx,thini,params,algo,bounds,aguess)
% [thhat,covh,logli,thini,scl,params,eflag,oput,grd,hes,Hk,k,options,bounds]=...
%          MLEOSL(Hx,thini,params,algo,bounds,aguess)
%
% Performs a maximum-likelihood estimation for SINGLE FIELDS as in
% Olhede & Simons (2013) by minimization using FMINUNC/FMINCON.
%
% INPUT:
%
% Hx       Real matrix with the field that is to be analyzed [whichever]
% thini    An unscaled starting guess for the parameter vector with elements:
%          [          s2 nu rho] - see SIMULOSL. If you leave this
%          value blank, then you will work from the perturbed "aguess"
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
%
% OUTPUT:
%
% thhat    The maximum-likelihood estimate of the vector [scaled]:
%          [s2 nu rho], in Nm, and "nothing", see SIMULOSL
% covh     The unscaled covariance from the (blurred) numerical Hessian at the estimate
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
% At least 'demo1' has been tested to run in an SPMD loop! Files are
% opened in append mode - except for thzro, which will only reflect one lab.
%
% EXAMPLE:
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
% Last modified by fjsimons-at-alum.mit.edu, 09/18/2016

% NEED TO CHANGE THE k(~~k) to proper accounting for kiso

if ~isstr(Hx)
  defval('algo','con')
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
                     [var(Hx)*100 8.00  max(5e5,min(dydx.*NyNx))],... % Upper bounds
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
  disp(sprintf(sprintf('\n%s : %s\n',str0,repmat(str2,size(thini))),...
	       'Starting theta',thini))

  % If you brought in your own initial guess, need an appropriate new scale
  if ~isempty(inputname(2)) || any(aguess~=thini)
    scl=10.^round(log10(abs(thini)));
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,size(scl))),...
		 'Scaling',scl))
  end
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

  % Modify to demean
  disp('DEMEAN BOTH DATA SETS')
  Hx(:,1)=Hx(:,1)-mean(Hx(:,1));
  % Let us NOT demean and see where we end up...

  % Turn the observation vector to the spectral domain
  % Watch the 2pi in SIMULOSL
  Hk(:,1)=tospec(Tx(:).*Hx(:,1),params)/(2*pi);

  NN=200;
  % And now get going with the likelihood using Hk
  % [ off|iter|iter-detailed|notify|notify-detailed|final|final-detailed ] 
  % Should probably make the tolerances relative to the number of k points
  options=optimset('GradObj','off','Display','on',...
		   'TolFun',1e-11,'TolX',1e-11,'MaxIter',NN,...
		   'LargeScale','off');
  % The 'LargeScale' option goes straight for the line search when the
  % gradient is NOT being supplied.

  % Set the parallel option to (never) use it for the actual optimization
  % Doesn't seem to do much when we supply our own gradient
  options.UseParallel='always';

  if blurs<2
    % Use the analytical gradient in the optimization, rarely a good idea
    % options.GradObj='on';
    if xver==1
      % Definitely good to check this once in a while
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
      % Could here compare to our own estimates of grad and hes!
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
    if xver==1
      disp(sprintf('%8.3gs per %i iterations or %8.3gs per %i function counts',...
		   ts/oput.iterations*100,100,ts/oput.funcCount*1000,1000))
    else
      disp(sprintf('\n'))
    end
  catch
    % If something went wrong, exit gracefully
    varargout=cellnan(nargout,1,1);
    return
  end

  % Only if extra verification is invoked, and no blurring was done,
  % since then the 'hes' that comes out of FMINUNC or FMINCON will be
  % based on unblurred data, and comparisons with analytics make sense
  if xver==1 && abs(blurs)<2
    % We have an analytic unblurred Hessian, so check the scaled versions
    H=Hessiosl(k,thhat.*scl,Hk).*[scl(:)*scl(:)'];
    % Calculate the theoretical covariance and Fisher matrices 
    [covF,FF]=covthosl(thhat,k,scl); 
    
    % So compare the analytic Hessian with the numerical Hessian and with
    % the Hessian expectation, which is the Fisher, all evaluated at the
    % estimate, and scaled. Note that we use the NEGATIVE LOG-LIKELIHOOD
    disp(' ')
    disp('Analytic Hessian       Observed -Hessian   Analytic Fisher')
    [trilos(H) trilos(hes) trilos(FF)]
    
    % This is the covariance as predicted from the analytic Hessian
    covH=hes2cov(H,scl,length(k(~~k))/2);
    
    % This is the covariance as calculated from the numerical Hessian
    covh=hes2cov(hes,scl,length(k(~~k))/2);
    disp('Analytic H-covariance  Observed H-covariance   Analytic F-covariance')
    [trilos(covH) trilos(covh) trilos(covF)]
    
    % how about the gradient? compare to grd./scl'
    grobs=-nanmean(gammakosl(k,thhat.*scl,params,Hk))';
    disp('Gradient ratios')
    [grobs./(grd./scl')]    
  end

  % This is the entire-plane estimate (hence the factor 2!) The
  % covariance as calculated from the numerically observed blurred Hessian
  covh=hes2cov(-hes,scl,length(k(~~k))/2);

  % Talk!
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(sprintf(sprintf('%s : %s\n ',str0,repmat(str2,size(thhat))),...
	       'Obs-Hessian std',sqrt(diag(covh))))

  % Generate output as needed
  varns={thhat,covh,logli,thini,scl,params,eflag,oput,grd,hes,Hk,k,options,bounds};
  varargout=varns(1:nargout);
elseif strcmp(Hx,'demo1')
  % If you run this again on the same date, we'll just add to THINI and
  % THHAT but you will start with a blank THZERO. See 'demo2'
  % How many simulations? The SECOND argument after the demo id
  defval('thini',[]);
  N=thini; clear thini
  defval('N',500)
  more off
  % What th-parameter set? The THIRD argument after the demo id
  defval('params',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=params; clear params
  % What fixed-parameter set? The FOURTH argument after the demo id
  defval('algo',[]);
  % If there is no preference, then that's OK, it gets taken care of
  params=algo; clear algo
    
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
    [thhat,covh,logli,thini,scl,p,e,o,gr,hs,Hk,k,ops,bnds]=...
	mleosl(Hx,[],p,[],[],th0);
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
		var(Hx),momx,covh)
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

    % This is the average of the Hessians, should be close to the Fisher
    avH=avH.*[sclth0(:)*sclth0(:)']/good;

    % Now compute the theoretical covariance and scaled Fisher
    [covF,F]=covthosl(th0./sclth0,k,sclth0);
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
  [th0,thhats,p,covF,covHav,covHts,~,~,~,~,momx]=osload(datum);

  % Report the findings of all of the moment parameters
  disp(sprintf('m(m(Xk)) %f m(v(Xk)) %f m(magic) %s v(magic) %f',...
	      mean(momx),var(momx(:,end))))
  
  % Plot it all - perhaps some outlier selection?
  [ah,ha]=mleplos(trimit(thhats,99,1),th0,covF,covHav,covHts,[],[],p,...
                  sprintf('MLEOSL-%s',datum));
  % Should report this below the number of data in the title

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
  
  % Take a look a the distribution of the momx?
  % With our definitions, the variance predicted (RB X p 51) should be
  % doh, that is the skewness of a chi-squared - just sayin'.
  % We don't change the number of degrees of freedom! If you have used
  % twice the number, and given half correlated variables, you do not
  % change the variance, that is the whole point. Unlike in COVTHOSL
  % where you make an analytical prediction that does depend on the
  % number and which therefore you need to adjust.
  k=knums(p); varpred=8/[length(k(~~k))];
  figure
  subplot(211)
  histfit(momx(:,3)); [m,s]=normfit(momx(:,3));
  subplot(212)
  qqplot((momx(:,3)-1)/sqrt(varpred)); axis image; grid on
  refline(1,0)
  % Could also do, as these quantities should be very close of course
  % qqplot(momx(:,2),momx(:,3)); axis image; refline(1,0); grid on
  % Then use NORMTEST to ascertain the veracity... don't bother with the
  % Nyquist wavenumbers, there will be very few, but take out the zero
  % Predicted expected value is one.
  [a,b,c,d]=normtest(momx(:,3),1,varpred);
  disp(sprintf('%i%% rejected at the %i%% confidence level',...
               round(c),round(d*100)))  
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
