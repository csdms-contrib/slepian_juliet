function [th0,thhats,params,covF,covHav,covHts,E,v,obscov,sclcovF,momx]=osload(datum,perc)
% [th0,thhats,params,covF,covHav,covHts,E,v,obscov,sclcovF,momx]=OSLOAD(datum,perc)
%
% Loads all four of the output diagnostic files produced by the suite of
% programs following Simons & Olhede (2013). 
%
% INPUT:
%
% datum        The file identifier
% perc         The percentage preserved by calling TRIMIT to the estimated parameters
%
% OUTPUT:
%
% th0          The true parameter vector
% thhats       The estimates from the inversion
% params       The parameter structure of the simulation
% covF         The covariance matrix based on the Fisher matrix at the truth
% covHav       The covariance matrix based on the average Hessian matrix
% covHts       The covariance matrix based on one randomly picked Hessian matrix
% E            Young's modulus - for conversion to Te only
% v            Poisson's ratio - for conversion to Te only
% obscov       Scaled observed covariance (ones on the diagonal)
% sclcovF      Scaled theoretical covariance (ones on the diagonal)
% momx         Moment parameters of the simulation results
%
% SEE ALSO:
%
% OSOPEN, DIAGNOS, TRIMIT, MLEOS etc
%
% Last modified by fjsimons-at-alum.mit.edu, 06/27/2016

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

% Then load and display what we have
fid=fopen(f1,'r');
[th0,params,sclth0,avhsz,Fisher,covF,nh]=osrzero(fid,np);
fclose(fid);

% Load the optimization diagnostics, which should duplicate f2 and f3
[thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covh]=diagnos(f4,pwd,np); 

try
  % Could be off; know that DIAGNOS is the file of record
  difer(thhat(:,1:np)-thhats,[],[],NaN)
  difer(thini        -thinis,[],[],NaN)
catch
    warning('THINI, THHAT and DIAGN do not agree (DIAGN will be used)')
end

% Here you might confirm that the hokey variances are on average further
% from the truth, negatively biased, than the mle estimates. 

% Bring all of them on a common scaling
if sum(sum(abs(diff(scl,1)),1)) || sum(abs(sclth0-unique(scl,'rows')))
  % Remember the way the scaling works for Hessians (which are inverse th0)
  for index=1:size(hes,1)
    hes(index,:)=trilos(trilosi(...
        hes(index,:))./[scl(index,:)'*scl(index,:)]...
                        .*[sclth0(:)*sclth0(:)']);
  end
end

% So THIS actually observed average is to be compared with the theory
% Note that this is sometimes very strongly affected by outliers;
% shouldn't give it too much credence. Try the median?
avhs=median(hes,1);

% This is the covariance derived from the grand average over all runs
% collected in the DIAGN file, expressed with a common scaling
k=knums(params);
covHav=hes2cov(trilosi(avhs),sclth0,length(k(~~k))*2);
% This is the "partial average" over the last nh runs as reported in the
% OSWZERO file, expressed with the same scaling again... wasn't the median
% If it was a "blank" filler shot at the end, ignore warnings
warning off MATLAB:singularMatrix
covHavz=hes2cov(trilosi(avhsz),sclth0,length(k(~~k))*2);
warning on MATLAB:singularMatrix

% Also pick some random ones from the list of Hessian-derived covariances
pix=randi(size(covh,1));
covHts=trilosi(covh(pix,:));
% Tell us what the three picked values were!
disp(sprintf('\n%s picking %i %g %i for the blue display\n',...
             upper(mfilename),thhats(pix,1:3)))
% The theoretical AT the random one isn't the same as the observed one at
% the random one, and THAT is what the blue line is.
covthpix=covthosl(thhats(pix,:),k,1);
% This should be the same as covF
covth0=covthosl(th0,k,1);

% Take a good look at the second one, just for kicks
disp(sprintf('black nu %f',indeks(sqrt(diag(covth0)),2)))
disp(sprintf('blue  nu %f',indeks(sqrt(diag(covthpix)),2)))

% fix the theoretical covariance to have physical scaling AT THE FIRST WRITING

% Remember the avhs is actually coming from the diagnostic file; below is untrimmed
osdisp(th0,thhat(:,1:np),size(hes,1),avhs,Fisher,covF)

if np>3
  % For display and later output only
  omukm=round(mean(DtoTe(thhat(:,1),E,v)/1000));
  ostkm=round(std(DtoTe(thhat(:,1),E,v)/1000));
  disp(sprintf('Average value of all  the Te estimates:  %4i km',omukm))
  disp(sprintf('Standard deviation of the Te estimates:  %3.2f km',ostkm))
end

% Easiest to compare the scaled versions
% Rescale the covariance matrix
%sclcovF=    covF./[sclth0(:)*sclth0(:)'];
%obscov=cov(thhat)./[sclth0(:)*sclth0(:)'];

% Actually, best to scale so the diagonal has unit variance
sclcovF=covF./[diag(sqrt(covF))*diag(sqrt(covF))'];

% How about we plot the observed covariance matrix instead
obscov=cov(thhat(:,1:np));
obscov=obscov./[diag(sqrt(obscov))*diag(sqrt(obscov))'];
