function varargout=sgp(params,Cmn)
% varargout=SGP(params,Cmn)
%
% Simulating stationary Gaussian field over a two-dimensional grid
%
% INPUT:
%
% params  parameters to the data grid, i.e. at least
%      params.dydx with the grid spacings
%      params.NyNx with the grid dimensions
% Cmn    scalar function handle to the covariance function
%         with one two-vector input C([m n])
% 
% OUTPUT:
%
%  f1,f2   two statistically independent fields over the mxn grid;
%  tx,ty   vectors so the field can be plotted via IMAGESC(tx,ty,f1)
%
% EXAMPLES:
%
%% Exponential covariance function, specified inline
% Cmn=@(h)([1-h(1)^2/th(1)^2-h(1)*h(2)/(th(2)*th(1))-h(2)^2/th(2)^2]...
%      *exp(-[h(1)^2/th(1)^2+h(2)^2/th(2)^2]));
%% Matern covariance, specified via 
% th0=[10 2.4 8]; params.NyNx=[512 384]; params.dydx=[1 1];
% Cmn=@(h) cov_matern_1(h,th,params.dydx);
%% Matern covariance, specified via MATERNOSY
% Cmn=@(h) maternosy(sqrt([h(1)*params.dydx(1)]^2+[h(2)*params.dydx(2)]^2),th0)
%% Perform the calculation
% sgp(params,Cmn)
%
%% Compare blurred simulation versus circulant-embedding simulation
% ah=krijetem(subnum(2,2));
% [Hx1,th0,params1,k,Hk1,Sb1]=simulosl; axes(ah(1)); imagesc(v2s(Hx1));
% title(sprintf('blurs %i',params1.blurs))
% params2=params1; params2.blurs=Inf;
% [Hx2,th0,params2,k,Hk2]=simulosl(th0,params2); axes(ah(2)); imagesc(v2s(Hx2));
% title(sprintf('blurs %i',params2.blurs))
%% This needs to absolutely agree
% axes(ah(3)); imagesc(v2s(log10(abs(Hk1).^2)));
%% axes(ah(3)); imagesc(v2s(log10(Sb1)));
% axes(ah(4)); imagesc(v2s(log10(abs(Hk2).^2)));
%
% Written by Arthur Guillaumin, 10/27/2017
% Last modified by fjsimons-at-alum.mit.edu, 08/16/2021

%% Reference:
% Kroese, D. P. & Botev, Z. I. (2015). Spatial Process Simulation.
% In Stochastic Geometry, Spatial Statistics and Random Fields
% DOI: 10.1007/978-3-319-10064-7_12 (pp. 369-404)

% M,N     row and column dimensions of the 2-D spatial output grid
%          (note that size of the covariance matrix is MNxMN)
M=params.NyNx(1);
N=params.NyNx(2);

% Create grid for spatial field
tx=[0:N-1]; 
ty=[0:M-1];
% Initialize auxiliary matrices
[Rows,Cols]=deal(zeros(M,N));
% Sample covariance function at grid points
for i=1:N 
  for j=1:M
    % Rows of blocks of cov matrix
    Rows(j,i)=Cmn([tx(i)-tx(1) ty(j)-ty(1)]);
    % Columns of blocks of cov matrix
    Cols(j,i)=Cmn([tx(1)-tx(i) ty(j)-ty(1)]);
  end
end
% FJS Really need to vectorize the above a bit, shouldn't I

% Create the first row of the block circulant matrix with circular blocks
% and store it as a matrix suitable for FFT2
BlkCirc_row=[Rows              Cols(:,end:-1:2);
             Cols(end:-1:2,:)  Rows(end:-1:2,end:-1:2)];
% Compute eigenvalues
lam=real(fft2(BlkCirc_row))/(2*M-1)/(2*N-1);
if abs(min(lam(lam(:)<0)))>1e-15
  error('Could not find positive definite embedding!')
else
  lam(lam(:)<0)=0; 
  lam=sqrt(lam);
end
% Generate field with covariance given by block circular matrix
F=fft2(lam.*complex(randn(2*M-1,2*N-1),randn(2*M-1,2*N-1)));
% Extract subblock with desired covariance
F=F(1:M,1:N); 
% Two independent fields with desired covariance
f1=real(F); 
f2=imag(F);

% Optional output
varns={f1,f2,tx,ty};
varargout=varns(1:nargout);

if nargout==0
  % Plot if no output requested
  imagesc(tx,ty,f1)
  colormap bone
  axis image
end

