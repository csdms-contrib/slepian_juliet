function [H,cH]=Hessiosl(k,th,params,Hk,xver)
% [H,cH]=HESSIOSL(k,th,params,Hk,xver)
%
% Calculates the entries in the Hessian matrix of Olhede & Simons (2013) for
% the Whittle-likelihood under the UNIVARIATE ISOTROPIC MATERN model, after
% wavenumber averaging. No consideration is given to the zero wavenumber,
% see LOGLIOSL.
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
% xver     Excessive verification [1 or 0]
%
% OUTPUT:
%
% H        The full-form Hessian matrix, a symmetric 3x3 matrix
% cF       The uniquely relevant elements listed in this order:
%          [1] Hs2s2   [2] Hnunu  [3] Hrhorho
%          [4] Hs2nu   [5] Hs2rho [6] Hnurho
%
% SEE ALSO:
%
% GAMMIOSL, FISHIOSL, HES2COV, TRILOS, TRILOSI
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=Fishiosl(k,th0); 
% G=gammiosl(k,th0,p,Hk);
% H=Hessiosl(k,th0,p,Hk);
% % On average, F and H should be close
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2016

% Early setup exactly as in FISHIOSL
defval('xver',1)

% Exclude the zero wavenumbers
Hk=Hk(~~k);
k=k(~~k);

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;
% The number of wavenumbers
lk=length(k(:));

% First compute the auxiliary parameters
[mth,mththp]=mAosl(k,th,xver);

% Extract the needed parameters of the simulation variables
blurs=params.blurs;
NyNx=params.NyNx;

% We need the power spectrum and its ratio to the observations 
% See LKOSL for the detailed explanations of these procedures
switch blurs
 case {0,1}
    S=maternos(k,th);
 otherwise
  if blurs>1
    S=bluros(maternos(knums(params,1),th),params);
  else
    S=blurosy(th,params);
  end
end

% The average of Xk needs to be close to one as will be tested 
Xk=abs(Hk).^2./S;

keyboard

% Allocate arrays for good measure; no cell since they're all k-dependent
Hk=nan(lk,npp);
Hf=nan(1,npp);
cH=nan(np,np);

% Hkthth
for j=1:3
  Hk(:,j)=-mththp{j}-[mth{j}.^2-mththp{j}].*Xk;
end

% All further combinations Hkththp
jcombo=nchoosek(1:np,2);
for j=1:length(jcombo)
  Hk(:,np+j)=-mththp{np+j}-[mth{jcombo(j,1)}.*mth{jcombo(j,2)}-mththp{np+j}].*Xk;
end

keyboard

% Now perform the averaging over all wavenumbers
Hf=nanmean(Hk,1);

% The full Hessian matrix
F=trilosi(cF);

