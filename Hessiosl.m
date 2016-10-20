function F=Hessiosl(k,th,params,Hk)
% F=HESSIOSL(k,th,params,Hk)
%
% Calculates the entries in the Fisher matrix of Olhede & Simons (2013) for
% the Whittle-likelihood estimate of the SINGLE-FIELD Matern model, after
% wavenumber averaging. When blurred, no consideration is given to the
% zero wavenumber, see LOGLIOSL. Zero-wavenumber excluded.
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
%
% OUTPUT:
%
% F        A full-form Hessian matrix [not scaled]
%
% SEE ALSO:
%
% FISHERKOSL, HES2COV
%
% EXAMPLE:
% 
% [~,th0,p,k,Hk]=simulosl([],[],1);
% F=Fishiosl(k,th0); H=Hessiosl(k,th0,p,Hk);
%% On average, these two should be close
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2016

% Early setup exactly as in FISHIOSL
defval('xver',1)

% Exclude the zero wavenumbers
k=k(~~k);

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;
% The number of wavenumbers
lk=length(k(:));

% First compute the "means" parameters, one per parametre
mth=mAosl(k,th,xver);

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

% Extract the parameters from the input
s2=th(1);
nu=th(2);
rho=th(3);

% Here be the derivat
% A variable that is also needed
vpiro=4/pi^2/rho^2;
avark=nu*vpiro+k(:).^2;
% Here is the matrix of derivatives of m, diagonals first, checked with 
% syms s2 nu rho k avark vpiro pi
% Checked dms2/ds2
mx{1}=-1/s2^2;
% Checked dmnu/dnu
mx{2}=2/nu-(nu+1)/nu^2-2*vpiro./avark+(nu+1)*vpiro^2./avark.^2;
% Checked dmrho/drho
mx{3}=2*nu*(1-3*(nu+1)*vpiro./avark+2*nu*(nu+1)*vpiro^2./avark.^2)/rho^2;
% Here the cross-terms, which are verifiably symmetric
mx{4}=0;
mx{5}=0;
% This I've done twice to check the symmetry, checked dmrhodnu or dmnudrho
mx{6}=2/rho*(-1+[2*nu+1]*vpiro./avark-nu*[nu+1]*vpiro^2./avark.^2);


% Allocate arrays for good measure; no cell since they're all k-dependent
Fk=nan(lk,npp);
Ff=nan(1,npp);
F =nan(np,np);

% Fkthth
for j=1:3
  Fk(:,j)=-mx{j}-[m{j}.^2-mx{j}].*Xk;
end

% All further combinations Fkththp
jcombo=nchoosek(1:np,2);
for j=1:length(jcombo)
  Fk(:,np+j)=-mx{np+j}-[m{jcombo(j,1)}.*m{jcombo(j,2)}-mx{np+j}].*Xk;
end

% Don't take out the zero wavenumber when there is no blurring, otherwise
% the score won't verify if you should choose to do so. On the whole, you
% won't want to run unblurred anything, so it won't be a big deal. 
if blurs~=0
  % Trouble is at the central wave numbers, we take those out 
  % Find the zero wavenumber
  kzero=sub2ind(NyNx,floor(NyNx(1)/2)+1,floor(NyNx(2)/2)+1);
  % Check that this is correctly done
  difer(k(kzero),[],[],NaN)
  % Behavior is rather different if this is NOT done... knowing that it
  % will not blow up but rather be some numerically large value
  Fk(kzero)=NaN;
end

% Now perform the averaging over all wavenumbers
Ff=nanmean(Fk,1);

% Ff     The 6-entry Hessian matrix, in this order:
%          [1] Fs2s2   [2] Fnunu  [3] Frhorho
%          [4] Fs2nu   [5] Fs2rho [6] Fnurho

% The upper half of full Hessian matrix
F(1,1)=Ff(1);
F(2,2)=Ff(2);
F(3,3)=Ff(3);

% These will be the covariances of D with others
F(1,2)=Ff(4);
F(1,3)=Ff(5);

% These will be the covariances of f2 with others
F(2,3)=Ff(6);

% Then symmetrize
F=[triu(F)'+triu(F)-diag(diag(F))];
