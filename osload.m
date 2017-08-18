function [th0,thhats,params,covFHh,covHav,thpix,E,v,obscov,sclcovF,momx,covthpix]=osload(datum,perc)
% [th0,thhats,params,covFHh,covHav,thpix,E,v,obscov,sclcovF,momx,covthpix]=OSLOAD(datum,perc)
%
% Loads ALL FOUR of the output diagnostic files produced by the suite of
% programs following Simons & Olhede (2013). Data files are, like,
% 'mleosl_thzro_16-Jun-2015-64-2', 'mleosl_thini_16-Jun-2015-64-2',
% 'mleosl_thhat_16-Jun-2015-64-2', 'mleosl_diagn_16-Jun-2015-64-2', etc.
%
% INPUT:
%
% datum        The file identifier, e.g. '16-Jun-2015-64-2'
% perc         The percentage preserved by calling TRIMIT to the estimated parameters
%
% OUTPUT:
%
% th0          The true parameter vector
% thhats       The maximum-likelihood estimates from the series of saved MLEOSL runs 
% params       The parameter structure of the SIMULOSL simulations
% covFHh       A covariance matrix for the estimate as written by OSWDIAG
% covHav       The covariance matrix based on the MEDIAN numerical Hessian matrix
% thpix        The example estimate, randomly picked up
% E            Young's modulus - for conversion to Te only
% v            Poisson's ratio - for conversion to Te only
% obscov       Scaled observed sample covariance (ones on the diagonal)
% sclcovF      Scaled Fisher covariance (ones on the diagonal)
% momx         Moment parameters of the simulation results
% covthpix     The covariance matrix based on the numerically derived
%              Hessian matrix at the randomly picked solution 
%
% SEE ALSO:
%
% OSOPEN, OSRDIAG, TRIMIT, MLEOS etc
%
% Last modified by fjsimons-at-alum.mit.edu, 08/18/2017

% Who called? Work this into the filenames
[~,n]=star69;

% Young's modulus 
defval('E',1.4e11);
% Poisson's ratio
defval('v',0.25);
% The percentile by which the parameters should be trimmed
defval('perc',100);

% You need to be in the proper directory to do this
f1=sprintf('%s_thzro_%s',n,datum);
f3=sprintf('%s_thhat_%s',n,datum);
f2=sprintf('%s_thini_%s',n,datum);
f4=sprintf('%s_diagn_%s',n,datum);

% Load the initial guesses
thinis=load(f2);
% Load the the estimates
thhats=load(f3);

% The number of parameters solved
np=size(thinis,2);

% Could be instructive to ascertain no patterns are in
% for i=1:np ; plot(thinis(:,i),thhats(:,i),'o'); pause; end

% Report what it is trying to readd
disp(sprintf('\n%s reading log files:\n\n%s\n%s\n%s\n%s',...
             upper(mfilename),f1,f2,f3,f4))

% Then load and display what we have
fid=fopen(f1,'r');
[th0,params,sclth0,avhsz,Fisher,covF,nh]=osrzero(fid,np);
fclose(fid);

% Load the optimization diagnostics, which should duplicate f2 and f3
[thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covFHh]=osrdiag(f4,pwd,np); 
try
  % Could be off; know that DIAGN is the file of record
  difer(thhat(:,1:np)-thhats,[],[],NaN)
  difer(thini        -thinis,[],[],NaN)
catch
    warning('THINI, THHAT and DIAGN do not agree (DIAGN will be used)')
end

% Here you might confirm that the hokey variances are on average further
% from the truth, negatively biased, than the mle estimates. 
matscl=[sclth0(:)*sclth0(:)'];

% Bring all of them on a common scaling
if sum(sum(abs(diff(scl,1)),1)) || sum(abs(sclth0-unique(scl,'rows')))
  % Remember the way the scaling works for Hessians (which are inverse th0)
  for index=1:size(hes,1)
    hes(index,:)=trilos(trilosi(...
        hes(index,:))./[scl(index,:)'*scl(index,:)].*matscl);
  end
end

% Below is the covariance derived from the "grand average" MEDIAN
% numerical Hessian over all runs collected in the DIAGN file, expressed
% with a common scaling. Note that this is sometimes very strongly
% affected by outliers; shouldn't give it too much credence.
avhs=median(hes,1);
k=knums(params);
df=length(k(~~k))/2; 
covHav=inv(trilosi(avhs)./matscl)/df;

% Below is the "partial average" MEAN over the last NH runs as reported in
% the OSWZERO file, expressed with the same scaling again...
% If it was a "blank" filler shot at the end, ignore, ignore warnings
covHavz=inv(trilosi(avhsz)./matscl)/df;

% A random pick from the set of maximum-likelihood estimates
pix=randi(size(thhats,1)); thpix=thhats(pix,1:3);
% Tell us what the three picked values were!
disp(sprintf('\n%s solution %i %g %i picked as an example\n',...
             upper(mfilename),thpix))

% A random pick from the set of numerical-Hessian derived covariances
covthpix=inv(trilosi(hes(pix,:))./matscl)/df;

% The below should be the same as covF since it came out of OSWZEROE, not at zero-k
[~,covth0]=fishiosl(k,th0,1);

% This needs to be close
diferm(covth0(:),covF(:),2)
% But the one taking out the wavenumber at zero is the right one
covF=covth0;

% Remember the avhs is actually coming from the diagnostic file; below is untrimmed
osdisp(th0,thhat(:,1:np),size(hes,1),avhs,Fisher,covHav)

if np>3
  % For display and later output only
  omukm=round(mean(DtoTe(thhat(:,1),E,v)/1000));
  ostkm=round(std(DtoTe(thhat(:,1),E,v)/1000));
  disp(sprintf('Average value of all  the Te estimates:  %4i km',omukm))
  disp(sprintf('Standard deviation of the Te estimates:  %3.2f km',ostkm))
end

% Actually, best to scale so the diagonal has unit variance
sclcovF=covF./[diag(sqrt(covF))*diag(sqrt(covF))'];

% How about we plot the observed sample covariance matrix instead
obscov=cov(thhat(:,1:np));
obscov=obscov./[diag(sqrt(obscov))*diag(sqrt(obscov))'];
