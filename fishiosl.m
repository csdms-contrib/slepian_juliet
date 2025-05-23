function [F,covF,cF]=fishiosl(k,th,params,xver)
% [F,covF,cF]=FISHIOSL(k,th,params,xver)
%
% Calculates the entries in the Fisher matrix of Olhede & Simons (2013, eq. 135)
% for the Whittle-likelihood under the UNIVARIATE ISOTROPIC MATERN model, after
% wavenumber averaging. No data are involved. Analytical expressions for the
% derivatives, see also LOGLIOSL. Zero-wavenumber excluded. No scaling asked or 
% applied. Exact blurring possible through MAOSL.

%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter [variance in unit^2]
%          th(2)=nu   The second Matern parameter [differentiability]
%          th(3)=rho  The third Matern parameter [range in m]
% params   Parameters of this experiment, the ones that are needed are:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 No wavenumber blurring
%                1 No wavenumber blurring, effectively
%                N Convolutional blurring, errors 
%               -1 Exact blurring
%              Inf Exact blurring, effectively
% xver     Excessive verification [0 or 1, which also computes F(k)]
%
% OUTPUT:
%
% F        The full-form Fisher matrix, a symmetric 3x3 matrix
% covF     A covariance estimate based on this matrix
% cF       The uniquely relevant Fisher elements listed in this order:
%          [1] F_s2s2   [2] F_nunu  [3] F_rhorho
%          [4] F_s2nu   [5] F_s2rho [6] F_nurho
%
% SEE ALSO: 
%
% GAMMIOSL, HESSIOSL, LOGLIOSL
% 
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=fishiosl(k,th0); g=gammiosl(k,th0,p,Hk); H=hessiosl(k,th0,p,Hk);
% round(abs((F+H)./F)*100) % should be single-digit small numbers
% [L,Lg,LH]=logliosl(k,th0,p,Hk);
% difer(Lg-g); difer(LH-H); % should be passing the test
%
% Last modified by fjsimons-at-alum.mit.edu, 04/8/2025
% Last modified by olwalbert-at-princeton.edu, 04/8/2025

% Early setup exactly as in HESSIOSL
defval('xver',1)

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% Exclude the zero wavenumbers
k=k(~~k);

% The number of nonzero wavenumbers
lk=length(k(:));

% First compute the auxiliary parameters
mth=mAosl(k,th,params,xver);

% Initialize
cF=nan(npp,1);

% Creative indexing - compare NCHOOSEK elsewhere
[i,j]=ind2sub([np np],trilos(reshape(1:np^2,np,np)));

% We're abusing the 'xver' switch to bypass saving wavenumber-dependencies
if xver==0
  % Attempt to do this otherwise for variance calculation, see also BLUROS
  % FJS [~,kz,W]=blurosy(th,params);
  % Do not save the wavenumber-dependent entities
  for ind=1:npp
    % Eq. (A60) in doi: 10.1093/gji/ggt056
    cF(ind)=mean(mth{i(ind)}.*mth{j(ind)});
    % Attempt to do this otherwise for variance calculation
    % FJS cFb(ind)=indeks(bsxfun(@times,mth{i(ind)},mth{j(ind)}'),':').*W(~~kz);
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

% Determine the degrees of freedom - could make sure to deduce this from
% the Hermiticity
df=lk/2;

% Construct the covariance matrix
covF=inv(F)/df;
