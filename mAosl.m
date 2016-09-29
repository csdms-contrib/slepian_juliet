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
% m        The "means" parameters which enter in the score, as a cell
% A        The "A" matrices which enter in the score, as a cell
%
% Last modified by fjsimons-at-alum.mit.edu, 06/22/2015

defval('xver',1)

% Extract the parameters from the input
s2=th(1);
nu=th(2);
rho=th(3);

% A variable that is also needed
avark=4*nu/pi^2/rho^2+k(:).^2;

% First calculate the "means" that depend on the wavenumbers
% Eq. (A25)
m{1}=1/s2;
% Eq. (A26)
m{2}=(nu+1)/nu+log(4*nu/pi^2/rho^2)...
     -4*(nu+1)/pi^2/rho^2./avark-log(avark);
% Eq. (A27)
m{3}=-2*nu/rho+8*nu/rho*(nu+1)/pi^2/rho^2./avark;

% Full output and extra verification etc
if nargout>1 || xver==1
  % Initialize
  A=cellnan(length(th),length(k(:)),3);
  % Obviously totally faking it here in this special case
  L=[ones(length(k(:)),1) zeros(length(k(:)),1) ones(length(k(:)),1)];

  A{1}=-repmat(m{1},length(k(:)),3);
  A{2}=-repmat(m{2},1,3);
  A{3}=-repmat(m{3},1,3);

  % Verification mode
  if xver==1
    % Check various identities at random wavenumbers 
    tracecheck(L,A,m,9,1)
  end
else
  A=NaN;
end

% Output
varns={m,A};
varargout=varns(1:nargout);
