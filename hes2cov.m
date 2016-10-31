function covH=hes2cov(H,scl,df)
% covH=HES2COV(H,scl,df)
%
% Turns a scaled Hessian into an unscaled covariance matrix according to
% asymptotic maximum-likelihood theory, see Simons & Olhede (2013). 
%
% INPUT:
% 
% H       A full-form Hessian matrix, e.g. from FMINUNC, FMINCON, or HESSIOSL
% scl     A vector with scaling factors applied to the parameter vector
% df      The number of degrees of freedom needed for final scaling, in
%         our case, this will be the number of independent wavenumbers
%         that went into the construction of the Whittle likelihood
%
% OUTPUT:
%
% covH    The covariance matrix based on this Hessian [not scaled]
%
% SEE ALSO:
%
% FISH2COV
%
% Last modified by fjsimons-at-alum.mit.edu, 10/13/2016

warning off MATLAB:nearlySingularMatrix

% Calculate this version of the covariance matrix
covH=inv(-H./[scl(:)*scl(:)'])/df;

warning on MATLAB:nearlySingularMatrix
