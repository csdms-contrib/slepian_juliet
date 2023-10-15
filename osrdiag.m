function [thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covX]=osrdiag(fname,ddir,np)
% [thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covX]=OSRDIAG(fname,ddir,np)
%
% Reads in a single file with diagnostics from MLEOS, MLEROS0, MLEROS.
%
% INPUT:
%
% ddir    A directory name, e.g. OLHEDE?/Simulations-lemaitre
% fname   A file name, e.g. 'diagn_21-Jun-2018'
% np      The number of parameters that are being solved for [default: 3]
%
% OUTPUT:
%
% See MLEOS, MLEROS0
%
% thhat    The maximum-likelihood estimate of the vector with elements:
%          [       s2 nu rho], in "nothing", see SIMULOSL
%          [D f2 r s2 nu rho], in Nm, and "nothing", see SIMULROS
%          [D f2   s2 nu rho], in Nm, and "nothing", see SIMULOS, SIMULROS0
% thini    The starting guess used in the optimization
% tseiter  Time taken, error code, number of iterations
% scl      The scaling applied as part of the optimization procedure
% L        The likelihood of this solution, by FMINUNC/FMINCON
% gam      The score, or first derivative of the likelihood, by FMINUNC/FMINCON
% hes      The Hessian, second derivative of the likelihood, by FMINUNC/FMINCON
% optis    The first-order optimality condition, by FMINUNC/FMINCON
% momx     The various moments of the quadratic piece of the likelihood
% covX     The covariance matrix estimate of the estimate
%
% SEE ALSO:
%
% OSOPEN, OSLOAD, OSWDIAG
%
% Last modified by fjsimons-at-alum.mit.edu, 10/14/2023

defval('ddir','/u/fjsimons/PROGRAMS/MFILES/olhede4')
defval('fname','mleosl_diagn_21-Jun-2018')

% The number of parameters to solve for
defval('np',3)
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% The number of sample variances that will be read - reset below
defval('nvar',2)

% Get a rough estimate of the number of estimates from the size
% You want this number to err on the side of being too large!
if np==3
  nsize=472; nvar=1;
elseif np==5
  nsize=560;
elseif np==6
  nsize=721;
end
ndim=ceil(fsize(fullfile(ddir,fname))/nsize);

% Initialize
tseiter=nan(ndim,3);
L=nan(ndim,1);
optis=nan(ndim,1);
momx=nan(ndim,3);
% Make room for the data sample variance(s)
thhat=deal(nan(ndim,np+nvar));
[thini,gam]=deal(nan(ndim,np));
hes=nan(ndim,npp);
covX=nan(ndim,npp);

% Rarely, in SPMD mode does the file get written too quickly and does a
% confusion between labs happen - look into the file and fix easily
% Read the contents

fid=fopen(fullfile(ddir,fname),'r');
for index=1:ndim
  try
    % The initial guess
    thini(index,:)=fscanf(fid,'%e',np);
  catch
    % Quit if you're done
    index=index-1;
    break
  end
  % The estimates, and the scaled spatial sample variance(s)
  thhat(index,:)=fscanf(fid,'%e',np+nvar);
  % Three diagnostics (time taken, exitflag, number of iterations)
  tseiter(index,:)=fscanf(fid,'%i',3); 
  try
    % The likelihood
    L(index)=fscanf(fid,'%e',1); 
    % disp(sprintf('L= %7.4f',L(index)))
  catch
    % Tell if you're wrong
    L(index-1)
    sprintf('%12.8g',thhat(index-1,3))
    break
  end
  % The first-order optimality criterion
  optis(index)=fscanf(fid,'%e',1);
  % The two moments of Xk, and the magic parameter
  momx(index,:)=fscanf(fid,'%e',3);
  % The scalings
  scl(index,:)=fscanf(fid,'%e',np);
  % The scores
  gam(index,:)=fscanf(fid,'%e',np);
  % The scaled Hessian elements at the solution
  hes(index,:)=fscanf(fid,'%f',npp);
  % An unscaled covariance estimate
  covX(index,:)=fscanf(fid,'%f',npp);
end
fclose(fid);

% Finetune the preallocation over time
if size(thhat,1)~=index
  disp(sprintf('\n%s allocated length %i, received length %i',...
               upper(mfilename),size(thhat,1),index))
end

% Trim the possibly wrongly preallocated arrays
thhat=thhat(1:index,:);
thini=thini(1:index,:);
tseiter=tseiter(1:index,:);
scl=scl(1:index,:);
L=L(1:index,:);
gam=gam(1:index,:);
hes=hes(1:index,:);
covX=covX(1:index,:);
optis=optis(1:index,:);
momx=momx(1:index,:);

% Put out
varns={thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covX};
varargout=varns(1:nargout);
