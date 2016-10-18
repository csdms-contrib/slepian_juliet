function F=Hessiosl(k,th,params,Hk)
% F=HESSIOSL(k,th,params,Hk)
%
% Calculates the entries in the Hessian matrix of Olhede & Simons (2013) 
% for the SINGLE-FIELD Matern model, post wavenumber averaging. 
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter, aka sigma^2 
%          th(2)=nu   The second Matern parameter 
%          th(3)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
%          NOTE: It's not going to be a great derivative unless you
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
% Last modified by fjsimons-at-alum.mit.edu, 10/18/2016

defval('xver',1)

% Usually we remove the zero wavenumber from consideration
if ~isempty(k(~k))
  disp(sprintf('%s zero wavenumber detected',upper(mfilename))); 
end

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% First compute the "means" parameters
m=mAosl(k,th,xver);

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

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
Fk=nan(length(k(:)),npp);
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

% Now perform the averaging over all wavenumbers
Ff=mean(Fk,1);

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

