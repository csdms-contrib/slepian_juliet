function [F,cF]=hessiosl(k,th,params,Hk,xver)
% [F,cF]=HESSIOSL(k,th,params,Hk,xver)
%
% Calculates the entries in the Hessian matrix of Olhede & Simons (2013) for
% the Whittle-likelihood under the UNIVARIATE ISOTROPIC MATERN model, after
% wavenumber averaging. Blurring is only approximately possible here, we
% work with analytical expressions for some of the derivatives, see
% LOGLIOSL. Zero-wavenumber excluded. No scaling asked or applied.
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter [variance]
%          th(2)=nu   The second Matern parameter [differentiability]
%          th(3)=rho  The third Matern parameter [range]
% params   A structure with AT LEAST these constants that are known:
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
%          NOTE: It's not going to be a great derivative unless you could
%          change MAOSL also. Still, the order of magnitude will be OK.
% Hk       A complex matrix of Fourier-domain observations
% xver     Excessive verification [0 or 1]
%
% OUTPUT:
%
% F        The full-form Hessian matrix, a symmetric 3x3 matrix
% cF       The uniquely relevant elements listed in this order:
%          [1] Fs2s2   [2] Fnunu  [3] Frhorho
%          [4] Fs2nu   [5] Fs2rho [6] Fnurho
%
% SEE ALSO:
%
% GAMMIOSL, FISHIOSL, HES2COV, TRILOS, TRILOSI
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); G=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be small numbers
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2016

% Early setup exactly as in FISHIOSL
defval('xver',1)

% Exclude the zero wavenumbers
Hk=Hk(~~k);
k=k(~~k);

% The number of parameters to solve for
np=length(th);
% The number of wavenumbers
lk=length(k(:));
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% First compute the auxiliary parameters
[mth,mththp]=mAosl(k,th,xver);

% Extract the needed parameters of the simulation variables
blurs=params.blurs;
NyNx=params.NyNx;

% We need the (blurred) power spectrum and its ratio to the observations
S=maternosp(k,th,params);

% The average of Xk needs to be close to one as will be tested 
Xk=abs(Hk).^2./S;

% Initialize
cF=nan(npp,1);

% Creative indexing - compare NCHOOSEK elsewhere
[i,j]=ind2sub([np np],trilos(reshape(1:np^2,np,np)));

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Do it all at once, don't save the wavenumber-dependent entities
  for ind=1:npp
    % Eq. (135) in doi: 10.1093/gji/ggt056
    cF(ind)=mean(-mththp{ind}-[mth{i(ind)}.*mth{j(ind)}-mththp{ind}].*Xk);
  end
elseif xver==1
  % Initialize; no cell since all of them depend on the wave vectors
  cFk=nan(lk,npp);
  % Do save the wavenumber-dependent entities
  for ind=1:npp
    cFk(:,ind)=-mththp{ind}-[mth{i(ind)}.*mth{j(ind)}-mththp{ind}].*Xk;
    % Eq. (135) in doi: 10.1093/gji/ggt056
    cF(ind)=mean(cFk(:,ind));
  end
end

% The full-form matrix
F=trilosi(cF);

% Can do an additional tracecheck here, using Option 2
