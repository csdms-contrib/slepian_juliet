function [L,gam,momx,vr]=logliosl(k,th,params,Hk,scl,xver)
% [L,gam,momx,vr]=LOGLIOSL(k,th,params,Hk,scl,xver)
%
% Calculates the full negative logarithmic likelihood and its
% derivatives, i.e. minus LKOSL and minus GAMMAKOSL averaged over
% wavenumber space. This is the function that we need to MINIMIZE! When
% blurred, no consideration is given to the zero wavenumber, see LKOSL. 
%
% INPUT:
%
% k        The wavenumbers at which these are being evaluated [1/m]
% th       The SCALED three-parameter vector argument
%          th(1)=s2   The first Matern parameter [variance in unit^2]
%          th(2)=nu   The second Matern parameter [differentiability]
%          th(3)=rho  The third Matern parameter [range in m]
% params   A structure with AT LEAST these constants that are known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%                -1 Blur likelihood using the exact BLUROSY procedure
%          kiso   wavenumber beyond which we are not considering the likelihood
% Hk       A [prod(params.NyNx)*1]-column of complex Fourier-domain observations
% scl      The SCALING factors applied, so that [scl.*thhat] is in units 
% xver     Excessive verification [0 or 1]
%
% OUTPUT:
%
% L        The loglihood, Lk averaged over all relevant wavenumbers
% gam      The score, averaged over all wavenumbers
% momx     Moments of the quadratic piece Xk over relevant wavenumbers;
%          the last one we use for a test on this being chi-squared
% vr       Variance of momx(3) under the null hypothesis
%
% SEE ALSO:
%
% HESSIOSL, FISHIOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 11/02/2016

defval('xver',1)

% Default scaling is none
defval('scl',ones(size(th)))

% Scale up the parameter vector for the proper likelihood and score
th=th.*scl;

% Here I build the protection that the three Matern parameters should be
% positive. I mirror them up! Thereby messing with the iteration path, but
% hey. It means we can use FMINUNC also.
th([1 2 3])=abs(th([1 2 3]));

% Filter, perhaps if-loop if extracting the Xk proves slower
[Lk,Xk]=Lkosl(k,th,params,Hk,xver);

if any(~isnan(params.kiso))
  Lk(k>params.kiso)=NaN;
  Xk(k>params.kiso)=NaN;
end

% Note: should we restrict this to the upper halfplane? or will mean do
% Get the likelihood at the individual wavenumbers; average
L=-nanmean(Lk);
% Attempt to reset if for some reason the whole thing failed
if isnan(L)
  L=1e100;
end

if nargout>=3
  % Extract the moments we'll be needing for evaluation later
  df=2;
  % First should be close to df/2, second close to df/2, third is like
  % the second except formally from a distribution that should be normal
  % with mean df/2 and variance 8/K --- the "magic" parameter
  momx=[nanmean(Xk) nanvar(Xk) nanmean([Xk-df/2].^2)];
end
if nargout>3
  % Compute the variance
  vr=8/sum(~isnan(Xk));
end

% I say, time to extract HESSIOSL here also?

% Get the appropriately scaled scores here
gam=gammiosl(k,th,params,Hk,xver).*scl(:);

% Print the trajectory, seems like one element at a time gets changed
%disp(sprintf('Current theta: %8.3g %8.3g %8.3g %8.3g %8.3g %8.3g',th))
