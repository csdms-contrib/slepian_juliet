function [L,g,H,momx,vr,Lk]=logliosl(k,th,scl,params,Hk,xver)
% [L,g,H,momx,vr,Lk]=LOGLIOSL(k,th,scl,params,Hk,xver)
%
% Calculates the full negative logarithmic likelihood and its derivatives as
% averaged over wavenumber space. This is the function that we need to
% MINIMIZE! No consideration is given to the zero wavenumbers.
%
% INPUT:
%
% k        The wavenumbers at which these are being evaluated [1/m]
% th       The three-parameter vector argument [scaled]
%          th(1)=s2   The first Matern parameter [variance in unit^2]
%          th(2)=nu   The second Matern parameter [differentiability]
%          th(3)=rho  The third Matern parameter [range in m]
% scl      The scaling factors applied, so that [scl.*thhat] is in units 
% params   A structure with AT LEAST these constants that are known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%                -1 Blur likelihood using the exact BLUROSY procedure
%          kiso   wavenumber beyond which we are not considering the likelihood
% Hk       A [prod(params.NyNx)*1]-column of complex Fourier-domain observations
% xver     Excessive verification [0 or 1, which also computes L(k)]
%
% OUTPUT:
%
% L        The logarithmic loglihood averaged over all nonzero wavenumbers
% g        The score averaged over all nonzero wavenumbers
% H        The Hessian averaged over all nonzero wavenumbers
% momx     Moments of the quadratic piece Xk over relevant wavenumbers;
%          the last one we use for a test on this being chi-squared
% vr       Variance of momx(3) under the null hypothesis
% Lk       The logarithmic loglihood prior to averaging over the wavenumbers
%
% SEE ALSO:
%
% GAMMIOSL, HESSIOSL, FISHIOSL
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); g=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be small numbers
% [L,Lg,LH]=logliosl(k,th0,1,p,Hk);
% difer(Lg-g); difer(LH-H); % should be passing the test
%
% Last modified by fjsimons-at-alum.mit.edu, 06/20/2018

defval('xver',1)

% Default scaling is none
defval('scl',ones(size(th)))

% Scale up the parameter vector for the proper likelihood and score
th=th.*scl;

% Here I build the protection that the three Matern parameters should be
% positive. I mirror them up! Thereby messing with the iteration path, but
% hey. It means we can use FMINUNC also.
th([1 2 3])=abs(th([1 2 3]));

% We need the (blurred) power spectrum and its ratio to the observations
[S,kk]=maternosp(th,params,xver);

% Exclude the zero wavenumbers
Hk=Hk(~~k);
S = S(~~kk); 
k = k(~~k);

% The statistics of Xk will be computed as momx, below, and tested in MLEOSL
Xk=hformos(S,Hk,[],xver);

% Radial wavenumber restrictions, these need to carry over elsewhere
if any(~isnan(params.kiso))
  ksel=k>params.kiso;
  disp('Need to still build this selection in...')
  disp('...depends if you want to look at Lk with the NaNs in there')
end

% Should make sure that this is real! Might remove REALIZE later
% Note: should we restrict this to the upper halfplane? or will MEAN do
% Should we use NANMEAN anywhere here? 

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Do it all at once, don't save the wavenumber-dependent entities
  % Eq. (A52) in doi: 10.1093/gji/ggt056
  L=-mean(-log(S)-Xk);
  Lk=NaN;
elseif xver==1
  % Do save the wavenumber-dependent entities
  Lk=realize(-log(S)-Xk);
  % Eq. (A52) in doi: 10.1093/gji/ggt056
  L=-mean(Lk);
end

% Attempt to reset if for some reason the whole thing failed
if isnan(L)
  % This may no longer be necessary% 
  L=1e100;
end

if nargout>=2
  % Get the appropriately scaled scores here
  g=gammiosl(k,th,params,Hk,xver).*scl(:);
end

if nargout>=3
  % Get the appropriately scaled Hessian values here
  H=hessiosl(k,th,params,Hk,xver).*[scl(:)*scl(:)'];
end

if nargout>=4
  % Extract the moments we'll be needing for evaluation later
  df=2;
  % First should be close to df/2, second close to df/2, third is like
  % the second except formally from a distribution that should be normal
  % with mean df/2 and variance 8/K --- the "magic" parameter
  momx=[nanmean(Xk) nanvar(Xk) nanmean([Xk-df/2].^2)];
end
if nargout>=5
  % Compute the variance
  vr=8/sum(~isnan(Xk));
end

% Print the trajectory, seems like one element at a time gets changed
%disp(sprintf('Current theta: %8.3g %8.3g %8.3g %8.3g %8.3g %8.3g',th))
