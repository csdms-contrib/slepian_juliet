function g=gammiosl(k,th,params,Hk,xver)
% g=GAMMIOSL(k,th,params,Hk,xver)
%
% Calculates the entries in the score matrix of Olhede & Simons (2013, eq. A. 53)
% for the Whittle-likelihood under the UNIVARIATE ISOTROPIC MATERN model,
% after wavenumber averaging. Blurring is only approximately possible here,
% we work with analytical expressions for some of the derivatives, see
% LOGLIOSL. Zero-wavenumber excluded. No scaling asked or applied.
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter [variance in unit^2]
%          th(2)=nu   The second Matern parameter [differentiability]
%          th(3)=rho  The third Matern parameter [range in m]
% params   A structure with AT LEAST these constants that are known:
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
%                Inf in which case it gets a hard reset to -1
%          NOTE: It's not going to be a great derivative unless you could
%          change MAOSL also. Still, the order of magnitude will be OK.
% Hk       A complex matrix of Fourier-domain observations
% xver     Excessive verification [0 or 1, which also computes g(k)]
%
% OUTPUT:
%
% g      The derivative of the log-likelihood, with elements
%          [1] g_s2   [2] g_nu  [3] g_rhor
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); g=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be small numbers
% [L,Lg,LH]=logliosl(k,th0,p,Hk);
% difer(Lg-g); difer(LH-H); % should be passing the test
% 
% Last modified by olwalbert-at-princeton.edu, 04/8/2024
% Last modified by fjsimons-at-alum.mit.edu, 04/8/2024

% params.blurs=Inf can only refer to spatial-domain generation and at
% this point we are already in the spectral domain; reset not returned
if isinf(params.blurs); params.blurs=-1; end

defval('xver',1)

% The number of parameters to solve for
np=length(th);

% We need the (blurred) power spectrum and its ratio to the observations
[S,kk]=maternosp(th,params,xver);

if xver==1
    diferm(kk(:),k(:))         
end

% Exclude the zero wavenumbers
Hk=Hk(~~k);
S = S(~~kk); 
k = k(~~k);

% The statistics of Xk are computed in LOGLIOSL
Xk=hformos(S,Hk,[],xver);

% First compute the auxiliary parameters
[mth,~,~,kkm]=mAosl(k,th,params,xver);

% The number of nonzero wavenumbers
lk=length(k(:));

% Initialize
g=nan(np,1);

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Do it all at once, don't save the wavenumber-dependent entities
  for ind=1:np
    % Eq. (A53) in doi: 10.1093/gji/ggt056
    g(ind)=-nanmean(-mth{ind}.*[1-Xk]);
  end
elseif xver==1 
  % Initialize; no cell since all of them depend on the wave vectors
  gk=nan(lk,np);
  % Do save the wavenumber-dependent entities
  for ind=1:np
    gk(:,ind)=mth{ind}.*[1-Xk];
    % Eq. (A53) in doi: 10.1093/gji/ggt056
    g(ind)=nanmean(gk(:,ind));
    % A somewhat redundant alternative way of computing these things
    diferm(gk(:,ind),mth{ind}-hformos(S,Hk,mth{ind}),7);
    % Pick up the additional necessities for a third way...
    [~,~,A]=mAosl(k,th,xver);
    diferm(gk(:,ind),mth{ind}+hformos(S,Hk,A{ind}),7);
  end
end
