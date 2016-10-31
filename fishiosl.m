function [F,cF]=fishiosl(k,th,xver)
% [F,cF]=FISHIOSL(k,th,xver)
%
% Calculates the entries in the Fisher matrix of Olhede & Simons (2013) for
% the Whittle-likelihood under the UNIVARIATE ISOTROPIC MATERN model, after
% wavenumber averaging. No blurring is possible here, since no data are
% involved and we work with analytical expressions for the derivatives, see
% LOGLIOSL. Zero-wavenumber excluded. No scaling asked or applied.
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter [variance]
%          th(2)=nu   The second Matern parameter [differentiability]
%          th(3)=rho  The third Matern parameter [range]
% xver     Excessive verification [0 or 1]
%
% OUTPUT:
%
% F        The full-form Fisher matrix, a symmetric 3x3 matrix
% cF       The uniquely relevant elements listed in this order:
%          [1] Fs2s2   [2] Fnunu  [3] Frhorho
%          [4] Fs2nu   [5] Fs2rho [6] Fnurho
%
% SEE ALSO: 
%
% GAMMIOSL, HESSIOSL, FISH2COV, TRILOS, TRILOSI
% 
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); G=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be small numbers
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2016

% Early setup exactly as in HESSIOSL
defval('xver',0)

% Exclude the zero wavenumbers
k=k(~~k);

% The number of parameters to solve for
np=length(th);
% The number of wavenumbers
lk=length(k(:));
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% First compute the auxiliary parameters
mth=mAosl(k,th,xver);

% Initialize
cF=nan(npp,1);

% Creative indexing - compare NCHOOSEK elsewhere
[i,j]=ind2sub([np np],trilos(reshape(1:np^2,np,np)));

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Do not save the wavenumber-dependent entities
  for ind=1:npp
    % Eq. (A60) in doi: 10.1093/gji/ggt056
    cF(ind)=mean(mth{i(ind)}.*mth{j(ind)});
  end
elseif xver==1
  % Initialize; some of them depend on the wave vectors, some don't
  cFk=cellnan([npp 1],[1 repmat(lk,1,5)],repmat(1,1,6));
  % Do save the wavenumber-dependent entities
  for ind=1:npp
    cFk{ind}=mth{i(ind)}.*mth{j(ind)};
    % Eq. (A60) in doi: 10.1093/gji/ggt056
    cF(ind)=mean(cFk{ind});
  end
end

% The full-form matrix
F=trilosi(cF);
    
