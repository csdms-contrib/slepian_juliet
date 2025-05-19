function [th0,thhats,params,covX,covavhs,thpix,E,v,obscov,sclcovX,momx,covXpix,covF0]=osload(datum,perc,n)
% [th0,thhats,params,covX,covavhs,thpix,E,v,obscov,sclcovX,momx,covXpix,covF0]=OSLOAD(datum,perc,n)
%
% Loads ALL FOUR of the output diagnostic files produced by the suite of
% programs following Simons & Olhede (2013). Data files are, like,
% 'mleosl_thzro_DD-MON-YYYY', 'mleosl_thini_DD-MON-YYYY'
% 'mleosl_thhat_DD-MON-YYYY', 'mleosl_diagn_DD-MON-YYYY', etc
%
% INPUT:
%
% datum        The file identifier, e.g. '16-Jun-2015-64-2'
% perc         The percentage preserved by calling TRIMIT to the estimated parameters
% n            A string if you do not want to accept the default calling-routine name
%
% OUTPUT:
%
% th0          The true parameter vector
% thhats       The maximum-likelihood estimates from the series of saved MLEOSL runs 
% params       The parameter structure of the SIMULOSL simulations
% covX         A covariance matrix estimate of the estimate as written by OSWDIAG
% covavhs      The covariance matrix based on the MEDIAN numerical Hessian matrix
% thpix        The example estimate, randomly picked up
% E            Young's modulus - for conversion to Te only
% v            Poisson's ratio - for conversion to Te only
% obscov       Scaled observed sample covariance (ones on the diagonal)
% sclcovX      Scaled version of a random pick of covX (ones on the diagonal)
% momx         Moment parameters of the simulation results
% covXpix      The covariance matrix based on the numerically derived
%              Hessian matrix at the randomly picked solution  
% covF0        The covariance matrix based on the Fisher matrix at the truth
%
% SEE ALSO:
%
% OSOPEN, OSRDIAG, TRIMIT, MLEOS etc
%
% Last modified by fjsimons-at-alum.mit.edu, 04/15/2025
% Last modified by olwalbert-at-princeton.edu, 04/15/2025

% Who called? Work this into the filenames
[~,nn]=star69;
defval('n',nn)

% Young's modulus 
defval('E',1.4e11);
% Poisson's ratio
defval('v',0.25);
% The percentile by which the parameters should be trimmed
defval('perc',99);

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

% Report what it is trying to read
disp(sprintf('\n%s reading log files:\n\n%s\n%s\n%s\n%s',...
             upper(mfilename),f1,f2,f3,f4))

% Then load and display what we have
fid=fopen(f1,'r');
[th0,params,sclth0,avhsz,F0,covF0,nh]=osrzero(fid,np);
fclose(fid);

% Load the optimization diagnostics, which should duplicate f2 and f3
[thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covX]=osrdiag(f4,pwd,np); 

try
  % Could be off; know that DIAGN is the file of record
  difer(thhat(:,1:np)-thhats,[],[],NaN)
  difer(thini        -thinis,[],[],NaN)
catch
    warning('THINI, THHAT and DIAGN do not agree (DIAGN will be used)')
end

% Here you might confirm that the hokey variances are on average further
% from the truth, negatively biased, than the MLE estimates. 
matscl=[sclth0(:)*sclth0(:)'];

% Bring all of them on a common scaling, if it should have differed...
if sum(sum(abs(diff(scl,1)),1)) || sum(abs(sclth0-unique(scl,'rows')))
  % Remember the way the scaling works for Hessians (which are inverse th0)
  for index=1:size(hes,1)
    hes(index,:)=trilos(trilosi(...
        hes(index,:))./[scl(index,:)'*scl(index,:)].*matscl);
  end
end

% Below is the covariance derived from the "grand average" MEDIAN
% numerical Hessian over all runs collected in the DIAGN file, expressed
% with a common scaling. Note that the average is very strongly
% affected by outliers; hence we pick the median.
avhs=median(hes,1);
k=knums(params);
df=length(k(~~k))/2; 
covavhs=inv(trilosi(avhs)./matscl)/df;

% Below is the "partial average" MEAN over the last NH runs as reported in
% the OSWZERO file, expressed with the same scaling again...
% If it was a "blank" filler shot at the end, ignore, ignore warnings
% covavhsz=inv(trilosi(avhsz)./matscl)/df;

% Now I am trimming
thhats=trimit(thhats,perc,1);

% A random pick from the set of maximum-likelihood estimates
pix=randi(size(thhats,1)); thpix=thhats(pix,1:3);
% Tell us what the three picked values were!
disp(sprintf('\n%s solution %i %g %i picked as an example\n',...
             upper(mfilename),thpix))

% A random pick from the set of numerical-Hessian derived covariances
covXpix=inv(trilosi(hes(pix,:))./matscl)/df;

% This should now be the same, the way the latest files have been written
pp=trilosi(covX(pix,:));
diferm(covXpix(:)./matscl(:),pp(:)./matscl(:),3)

% The below should be the same as covF0 since it came out of OSWZEROE, not at zero-k
[~,covth0]=fishiosl(k,th0,params,1);

% This needs to be identical except in a few digits compared to the size
diferm(covth0(:)./matscl(:),covF0(:)./matscl(:),3)
% But the one taking out the wavenumber at zero is the right one
covF0=covth0;

% Remember the avhs is actually the MEDIAN coming from the diagn file and
% not the AVERAGE coming from the thzero file, and the below is untrimmed
osdisp(th0,thhats(:,1:np),size(hes,1),avhs,F0,covavhs)

if np>3
  % For display and later output only
  omukm=round(mean(DtoTe(thhat(:,1),E,v)/1000));
  ostkm=round(std(DtoTe(thhat(:,1),E,v)/1000));
  disp(sprintf('Average value of all  the Te estimates:  %4i km',omukm))
  disp(sprintf('Standard deviation of the Te estimates:  %3.2f km',ostkm))
end

% Actually, best to scale so the diagonal has unit variance
sclcovX=covXpix./[diag(sqrt(covXpix))*diag(sqrt(covXpix))'];

% How about we return the observed sample covariance matrix instead
obscov=cov(thhats(:,1:np));
obscov=obscov./[diag(sqrt(obscov))*diag(sqrt(obscov))'];
