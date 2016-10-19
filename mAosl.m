function varargout=mAosl(k,th,xver)
% [m,A]=mAosl(k,th,xver)
%
% Calculates some auxiliary numbers and matrices for Olhede & Simons
% in the SINGLE-FIELD case. With some obvious redundancy for A.
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with elements [not scaled]:
%          s2  The first Matern parameter, aka sigma^2
%          nu  The second Matern parameter
%          rho The third Matern parameter
%
% OUTPUT:
%
% m        The "means" parameters for GAMMIOSL, HESSIOSL, as a cell
% A        The "A" matrices for GAMMIOSL, HESSIOSL, as a cell
%
% Last modified by fjsimons-at-alum.mit.edu, 10/19/2016

% Extra verification?
defval('xver',1)

% Extract the parameters from the input
s2=th(1);
nu=th(2);
rho=th(3);

% How many wavenumbers?
lk=length(k(:));

% Auxillary variable
avark=4*nu/pi^2/rho^2+k(:).^2;

% The "means",which still depend on the wavenumbers, sometimes
% Eq. (A25) in doi: 10.1093/gji/ggt056
m{1}=1/s2;
% Eq. (A26) in doi: 10.1093/gji/ggt056
m{2}=(nu+1)/nu+log(4*nu/pi^2/rho^2)...
     -4*(nu+1)/pi^2/rho^2./avark-log(avark);
% Eq. (A27) in doi: 10.1093/gji/ggt056
m{3}=-2*nu/rho+8*nu/rho*(nu+1)/pi^2/rho^2./avark;

% Full output and extra verification etc
if nargout>1 || xver==1
  % Initialize
  A=cellnan(length(th),lk,3);
  % Eq. (A28) in doi: 10.1093/gji/ggt056
  A{1}=-repmat(m{1},lk,3);
  % Eq. (A29) in doi: 10.1093/gji/ggt056
  A{2}=-repmat(m{2},1,3);  
  % Eq. (A30) in doi: 10.1093/gji/ggt056
  A{3}=-repmat(m{3},1,3);

  % Verification mode
  if xver==1
    % In this univariate case, we have some rather simple forms
    L=[ones(lk,1) zeros(lk,1) ones(lk,1)];
    % Check various identities at random wavenumbers 
    tracecheck(L,A,m,9,1)
  end
else
  A=NaN;
end

% Output
varns={m,A};
varargout=varns(1:nargout);
