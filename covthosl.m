function covth=covthosl(th,params,covm,ifinv)
% covth=COVTHOSL(th,params,covg,ifinv)
%
% Estimates the covariance matrix of the parameter estimate according to Eq. 36
% of Guillaumin et al 2022 as the matrix product of the Fisher information 
% matrix with the covariance matrix of the score
%
% INPUT:
% 
% th        Matern spectral parameters, [s2 nu rh]
% params    Parameter structure of the grid; at least p.NyNx and p.dydx must be
%           provided explicitly
% covm      A method specification for the covariance calculation
%           1 sampling
%           2 dfmtx
%           3 diagonals
% ifinv     Indicate which parameters were inverted for and require covariance
%           calculation
%
% OUTPUT:
%
% covth     Covariance of the parameter estimates
%
% SEE ALSO:
%
% COVGAMMIOSL
%
% EXAMPLES:
%
% th=[1.41 0.75 6];
% p=[];p.NyNx=[65 92];p.dydx=[1 1];p.blurs=Inf;ifinv=[1 1 1];
% covg=covgammiosl(th,p,1,ifinv);
% covth=covthosl(th,p,covm,ifinv);
%
% covth=covthosl(th,p,1,ifinv);
%
% Last modified by fjsimons-at-alum.mit.edu, 05/15/2025
% Last modified by olwalbert-at-princeton.edu, 05/15/2025

if covm==0
    covth=zeros(sum(ifinv),sum(ifinv));
    return
end

% Let's stick to easier problems for now
defval('ifinv',[1 0 1])

% If we did not bring in a precomputed variance of the score
% calculate it now via the specified method
covm=covgammiosl(th,params,1,ifinv);

% Calculate the Fisher information matrix for the grid
k=knums(params);
F=fishiosl(k,th,params);

% If the data were tapered and the periodogram was blurred, I think we need to
% modify how F is normalized in FISHIOSL (i.e., taking the mean of mth products)
if isfield(params,'taper') & numel(params.taper)>1 & isinf(params.blurs)
  nt=sum(params.taper,"all");
  nt=nt/prod(params.NyNx); 
  % Should be less than or equal to 1; if not, check out the taper
  if nt>1;keyboard;end
  F=F*nt;
end

% Degrees of freedom calculated according to FISHIOSL
df=length(k(:))/2;

% Only calculate the covariance for parameters that we inverted for
F=matslice(F,ifinv);
covF=inv(F);
covth=covF*covm*covF;

