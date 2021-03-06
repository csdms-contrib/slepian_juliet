function [Lk,Xk]=Lkosl(k,th,params,Hk,xver)
% [Lk,Xk]=Lkosl(k,th,params,Hk,xver)
%
% Computes the likelihood function for the SINGLE-FIELD isotropic Matern
% model in Olhede & Simons (2013). Always sets the zero wavenumber-value
% of Lk to NaN.
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
%          The zero wavenumber is present but excluded from consideration
% th       The parameter vector with elements [unscaled]:
%          s2  The first Matern parameter, aka sigma^2
%          nu  The second Matern parameter
%          rho The third Matern parameter
% params   A structure with AT LEAST these constants that are known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
% Hk       A [prod(params.NyNx)*1]-column vector of complex Fourier-domain observations
% xver     Excessive verification [0 or 1]
%
% OUTPUT:
%
% Lk       A one-column likelihood vector with the wavenumbers unwrapped
% Xk       A quadratic piece of it that gets used in the analysis of residuals
%
% SEE ALSO: 
%
% LOGLIOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 11/15/2016

% Extra verification?
defval('xver',1)

% Extract the needed parameters of the simulation variables
NyNx=params.NyNx;

% We need the (blurred) power spectrum and its ratio to the observations
S=maternosp(th,params,xver);

% The average of Xk needs to be close to one as will be tested 
warning off MATLAB:log:logOfZero
Xk=hformos(S,Hk,[],xver);
warning on MATLAB:log:logOfZero
% Should make sure that this is real! Don't take any chances
Lk=realize(-log(S)-Xk);

% Always remove the zero wavenumber from consideration
% Find the zero wavenumber
kzero=sub2ind(NyNx,floor(NyNx(1)/2)+1,floor(NyNx(2)/2)+1);
% Check that this is correctly done
difer(k(kzero),[],[],NaN)
% Behavior is rather different if this is NOT done... knowing that it
% will not blow up but rather be some numerically large value
Lk(kzero)=NaN;

