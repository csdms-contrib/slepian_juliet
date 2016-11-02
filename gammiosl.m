function G=gammiosl(k,th,params,Hk,xver)
% G=GAMMIOSL(k,th,params,Hk,xver)
%
% Calculates the entries in the score matrix of Olhede & Simons (2013) for
% the Whittle-likelihood under the UNIVARIATE ISOTROPIC MATERN model, after
% wavenumber averaging. Blurring is only approximately possible here, we
% work with analytical expressions for some of the derivatives, see
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
%          NOTE: It's not going to be a great derivative unless you could
%          change MAOSL also. Still, the order of magnitude will be OK.
% Hk       A complex matrix of Fourier-domain observations
% xver     Excessive verification [0 or 1]
%
% OUTPUT:
%
% G      The derivative of the log-likelihood, with elements
%          [1] Gs2s2   [2] Gnunu  [3] Grhorho
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); G=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be small numbers
% 
% Last modified by fjsimons-at-alum.mit.edu, 11/2/2016

defval('xver',1)

% Exclude the zero wavenumbers
Hk=Hk(~~k);
k=k(~~k);

% The number of parameters to solve for
np=length(th);
% The number of wavenumbers
lk=length(k(:));

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

% Initialize the wavenumber-dependent score vector
Gk=nan(lk,np);

% First compute the auxiliary parameters
if xver==0
  mth=mAosl(k,th,xver);
elseif xver==1
  [mth,~,A]=mAosl(k,th,xver);
end

% We need the (blurred) power spectrum and its ratio to the observations
[S,kk]=maternosp(th,params,xver);
% Exclude the zero wavenumbers
S=S(~~kk); 

% The average of Xk needs to be close to one as will be tested 
Xk=abs(Hk).^2./S;

% Initialize
G=nan(np,1);

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Do it all at once, don't save the wavenumber-dependent entities
  for j=1:np
    % Eq. (A53) in doi: 10.1093/gji/ggt056
    G(j)=-mean(-mth{j}.*[1-Xk]);
  end
elseif xver==1
  % Initialize; no cell since all of them depend on the wave vectors
  Gk=nan(lk,np);
  % Do save the wavenumber-dependent entities
  for j=1:np
    % Eq. (A53) in doi: 10.1093/gji/ggt056
    Gk(:,j)=-mth{j}.*[1-Xk];
    % Check for numerical indistinguishability ->>> MOVE TO ELSEWHERE
    gkalt1(:,j)=-mth{j}-hformos(S,A{j},Hk);
    gkalt2(:,j)=-mth{j}+hformos(S,mth{j},Hk);
    diferm(gkalt1(:,j),Gk(:,j),7);
    diferm(gkalt2(:,j),Gk(:,j),7);
  end
  % Now the wavenumber averaging
  G=-mean(Gk)';
end
