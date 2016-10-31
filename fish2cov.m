function [covF,F]=fish2cov(Fut,scl,df)
% [covF,F]=FISH2COV(Fut,scl,df)
%
% Turns (the upper-triangular portion of) a Fisher matrix into a covariance
% matrix and a scaled full-form Fisher matrix for Simons & Olhede (2013)
%
% INPUT:
% 
% Fut      (Upper-triangular portion of a) Fisher matrix [not scaled]
% scl      A vector with scaling factors applied to the parameter vector
% df       The number of degrees of freedom needed for final scaling, in
%          our case, this will be the number of independent wavenumbers
%          that went into the construction of the Whittle likelihood
%
% OUTPUT:
%
% covF     The covariance matrix based on this Fisher matrix [not scaled]
% F        The full-form Fisher matrix [scaled]
%
% SEE ALSO:
%
% HES2COV, TRILOSI, COVTHOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 06/22/2015

warning off MATLAB:nearlySingularMatrix

% Calculate this version of the covariance matrix
covF=inv(F)/df;

% Full form and in the scaled dimensions of the problem
F=[scl(:)*scl(:)'].*F;

warning on MATLAB:nearlySingularMatrix
