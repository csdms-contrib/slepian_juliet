function [H,covH,cH]=hessiosl(k,th,params,Hk,xver)
% [H,covH,cH]=HESSIOSL(k,th,params,Hk,xver)
%
% Calculates the entries in the Hessian matrix of Olhede & Simons (2013, eq. 132)
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
% xver     Excessive verification [0 or 1, which also computes H(k)]
%
% OUTPUT:
%
% H        The full-form Hessian matrix, a symmetric 3x3 matrix
% covH     A covariance estimate based on this Hessian matrix
% cH       The uniquely relevant Hessian elements listed in this order:
%          [1] H_s2s2   [2] H_nunu  [3] H_rhorho
%          [4] H_s2nu   [5] H_s2rho [6] H_nurho
%
% SEE ALSO:
%
% GAMMIOSL, FISHIOSL, LOGLIOSL
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); g=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be small numbers
% [L,Lg,LH]=logliosl(k,th0,p,Hk);
% difer(Lg-g); difer(LH-H); % should be passing the test
%
% Last modified by olwalbert-at-princeton.edu, 12/17/2024
% Last modified by fjsimons-at-alum.mit.edu, 12/17/2024

% params.blurs=Inf can only refer to spatial-domain generation and at
% this point we are already in the spectral domain; reset not returned
if isinf(params.blurs); params.blurs=-1; end

% Early setup exactly as in FISHIOSL
defval('xver',1)

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% We need the (blurred) power spectrum and its ratio to the observations
[S,kk]=maternosp(th,params,xver);

% Exclude the zero wavenumbers
Hk=Hk(~~k);
S=S(~~kk);
k=k(~~k);

% The number of nonzero wavenumbers
lk=length(k(:));

% The statistics of Xk are computed in LOGLIOSL
Xk=hformos(S,Hk,[],xver);

% First compute the auxiliary parameters
[mth,mththp]=mAosl(k,th,xver);

% Initialize
cH=nan(npp,1);

% Creative indexing - compare NCHOOSEK elsewhere
[i,j]=ind2sub([np np],trilos(reshape(1:np^2,np,np)));

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Do it all at once, don't save the wavenumber-dependent entities
  for ind=1:npp
    % Eq. (135) in doi: 10.1093/gji/ggt056
    cH(ind)=nanmean(-mththp{ind}-[mth{i(ind)}.*mth{j(ind)}-mththp{ind}].*Xk);
  end
elseif xver==1
  % Initialize; no cell since all of them depend on the wave vectors
  cHk=nan(lk,npp);
  % Do save the wavenumber-dependent entities
  for ind=1:npp
    cHk(:,ind)=-mththp{ind}-[mth{i(ind)}.*mth{j(ind)}-mththp{ind}].*Xk;
    % Eq. (135) in doi: 10.1093/gji/ggt056
    cH(ind)=nanmean(cHk(:,ind));
    % Can do an additional TRACECHECK here, using Option 2
    % And other tests, using HFORMOS, as in GAMMIOSL that for later
  end
end

% The full-form matrix
H=trilosi(cH);

% Determine the degrees of freedom - could make sure to deduce this from
% the Hermiticity
df=lk/2;

% Construct the covariance matrix
covH=inv(-H)/df;
