function varargout=mAosl(k,th,xver)
% [mth,mththp,A]=mAosl(k,th,xver)
%
% Calculates auxiliary derivative quantities for Olhede & Simons
% in the SINGLE-FIELD case. With some obvious redundancy for A.
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with elements [not scaled]:
%          th(1)=s2  The first Matern parameter [variance in unit^2]
%          th(2)=nu  The second Matern parameter [differentiability]
%          th(3)=rh  The third Matern parameter [range in m]
% xver     Excessive verification [0, 1 or 2]
%
% OUTPUT:
%
% mth      The first-partials "means" for GAMMIOSL, HESSIOSL, as a cell
% mththp   The second-partials parameters for GAMMIOSL, HESSIOSL, as a cell
% A        The "A" matrices for GAMMIOSL, HESSIOSL, as a cell
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2023

% Extra verification?
defval('xver',1)

% Extract the parameters from the input
s2=th(1);
nu=th(2);
rh=th(3);

% Auxilliary variable
avark=4*nu/pi^2/rh^2+k(:).^2;

% The "means", which still depend on the wavenumbers, sometimes
% The first derivatives of the logarithmic spectral density
% Eq. (A25) in doi: 10.1093/gji/ggt056
mth{1}=1/s2;
% Eq. (A26) in doi: 10.1093/gji/ggt056
mth{2}=(nu+1)/nu+log(4*nu/pi^2/rh^2)...
       -4*(nu+1)/pi^2/rh^2./avark-log(avark);
% Eq. (A27) in doi: 10.1093/gji/ggt056
mth{3}=-2*nu/rh+8*nu/rh*(nu+1)/pi^2/rh^2./avark;

if nargout>1
  % The second derivatives of the logarithmic spectral density
  vpiro=4/pi^2/rh^2;
  avark=nu*vpiro+k(:).^2;
  % Here is the matrix of derivatives of m, diagonals first, checked with 
  % syms s2 nu rh k avark vpiro pi
  % Checked dms2/ds2
  mththp{1}=-1/s2^2;
  % Checked dmnu/dnu
  mththp{2}=2/nu-(nu+1)/nu^2-2*vpiro./avark+(nu+1)*vpiro^2./avark.^2;
  % Checked dmrh/drh
  mththp{3}=2*nu*(1-3*(nu+1)*vpiro./avark+2*nu*(nu+1)*vpiro^2./avark.^2)/rh^2;
  % Here the cross-terms, which are verifiably symmetric
  mththp{4}=0;
  mththp{5}=0;
  % This I've done twice to check the symmetry, checked dmrhdnu or dmnudrh
  mththp{6}=2/rh*(-1+[2*nu+1]*vpiro./avark-nu*[nu+1]*vpiro^2./avark.^2);
else 
  mththp=NaN;
end

% Full output and extra verification
if nargout>2 || xver==1 || xver==2
  % How many wavenumbers?
  lk=length(k(:));
  
  % Initialize
  A=cellnan(length(th),lk,3);

  % Eq. (A28) in doi: 10.1093/gji/ggt056
  A{1}=-repmat(mth{1},lk,3);
  % Eq. (A29) in doi: 10.1093/gji/ggt056
  A{2}=-repmat(mth{2},1,3);  
  % Eq. (A30) in doi: 10.1093/gji/ggt056
  A{3}=-repmat(mth{3},1,3);

  % Verification mode
  if xver==1 || xver==2
    % In this univariate case, we have some rather simple forms
    % In the bivariate case, these are the nonzero lower-triangular
    % entries of a matrix that is another's Cholesky decomposition
    L=[ones(lk,1) zeros(lk,1) ones(lk,1)];
    tracecheck(L,A,mth,9,1)
    % How about the second TRACECHECK? We should be able to involve mththp
  end
else
  A=NaN;
end

% As much output as you desire
varns={mth,mththp,A};
varargout=varns(1:nargout);
