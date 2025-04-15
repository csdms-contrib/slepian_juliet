function covth=covthosl(th,params,covg,ifinv)
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
% covg      The variance of the score
% ifinv     Indicate which parameters were inverted for and require covariance
%           calculation
%
% OUTPUT:
%
% covth     Covariance of the parameter estimates; We are assuming a special
%           case for now so we only calculate a 2by2 for s2 and rh
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
% covth=covthosl(th,p,covg,ifinv);
%
% covth=covthosl(th,p,[],ifinv);
%
% Last modified by fjsimons-at-alum.mit.edu, 04/15/2025
% Last modified by olwalbert-at-princeton.edu, 04/15/2025

% We are calculating the variance of the estimates to incorporate wavenumber
% correlation when the periodogram is blurred. Confirm that we are using
% blurring
if ~isinf(params.blurs)
    params.blurs=Inf;
    warning('params.blurs was not, but now is, set to Inf')
end

% Let's stick to easier problems for now
defval('ifinv',[1 0 1])

% If we did not bring in a precomputed variance of the score
% calculate it now via sampling
if isempty(covg)
    covg=covgammiosl(th,params,1,ifinv);
end

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

% REMINDER: follow up on including the dof and normalization for taper
%
% Frederik normalizes by df in calculation of covF in FISHIOSL; Arthur does not
% in likelihood.DebiasedWhittle.fisher or 
% likelihood.DebiasedWhittle.variance_of_estimates
% I will leave out the degrees of freedom for now for comparison
% covFF=inv(F)/df;

covth=covF*jmat*covF;

% (sclh.*thhat,params,meth,ifinv);


% J=covgammiosl(sclh.*thhat,params,meth,ifinv);



% % obviously

% varianceoftheestimates

% % (mleoslcov could be )
