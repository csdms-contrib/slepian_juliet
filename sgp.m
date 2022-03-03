function varargout=sgp(params,Cnm)
% [f1,f2,x,y]=SGP(params,Cnm)
%
% Simulating stationary Gaussian field over a two-dimensional grid
%
% INPUT:
%
% params  parameters to the data grid, i.e. at least
%      params.dydx with the grid spacings
%      params.NyNx with the grid dimensions
% Cnm     VECTORIZED scalar function handle to the covariance function
%         with one two-vector input such that cov(x_t,y_s)=Cnm(t-s) is the
%         covariance function of a 2-dimensional stationary Gaussian field
% 
% OUTPUT:
%
%  f1,f2   two statistically independent fields over the mxn grid;
%  x,y   vectors so the field can be plotted via IMAGESC(x,y,f1)
%
% SEE ALSO:
%
% STATIONARY_GAUSSIAN_PROCESS2
%
% EXAMPLES:
%
% m=512; n=384;
% th0=[10 2.4 30]; params.NyNx=[512 384]; params.dydx=[10 10];
% rho2=@(h) maternosy(sqrt([h(:,1)*params.dydx(2)].^2+[h(:,2)*params.dydx(1)].^2),th0);
% figure(1)
% stationary_Gaussian_process2(m,n,rho2)
% t=title(sprintf('%s^2 %g, %s %g, %s %g','\sigma',th0(1),'\nu',th0(2),'\rho',th0(3))); movev(t,-m/30)
%
% figure(2)
% rho2=@(h) maternosy(sqrt([h(:,1)*params.dydx(2)].^2+[h(:,2)*params.dydx(1)].^2),th0);
% sgp(params,rho2)
% t=title(sprintf('%s^2 %g, %s %g, %s %g','\sigma',th0(1),'\nu',th0(2),'\rho',th0(3))); movev(t,-m/30)
%
%% Compare blurred simulation versus circulant-embedding simulation    
% figure(3)
% ah=krijetem(subnum(2,2));
% [Hx1,th0,params1,k,Hk1,Sb1]=simulosl; axes(ah(1)); imagesc(v2s(Hx1));
% title(sprintf('blurs %i',params1.blurs))
% params2=params1; params2.blurs=Inf;
% [Hx2,th0,params2,k,Hk2]=simulosl(th0,params2); axes(ah(2)); imagesc(v2s(Hx2));
% title(sprintf('blurs %i',params2.blurs))
%% These needs to absolutely agree
% axes(ah(3)); imagesc(v2s(log10(abs(Hk1).^2)));
%% axes(ah(3)); imagesc(v2s(log10(Sb1)));
% axes(ah(4)); imagesc(v2s(log10(abs(Hk2).^2)));
%
% Written by Arthur Guillaumin, 10/27/2017
% Last modified by fjsimons-at-alum.mit.edu, 02/28/2022

%% Reference:
% Kroese, D. P. & Botev, Z. I. (2015). Spatial Process Simulation.
% In Stochastic Geometry, Spatial Statistics and Random Fields
% DOI: 10.1007/978-3-319-10064-7_12 (pp. 369-404)

% M,N     row and column dimensions of the 2-D spatial output grid
%          (note that size of the covariance matrix is MNxMN)
M=params.NyNx(1);
N=params.NyNx(2);

% Create grid for spatial field in samples
tx=[0:N-1]; 
ty=[0:M-1];
% Initialize auxiliary matrices
[Rows,Cols]=deal(zeros(M,N));

% Sample covariance function at grid points
[TX,TY]=meshgrid(tx,ty);
% Remember that the function handle takes into account the sampling step
Rows=reshape(Cnm([ TX(:) TY(:)]),size(TX));
Cols=reshape(Cnm([-TX(:) TY(:)]),size(TX));

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

% Proper grid
x=tx*params.dydx(2);
y=ty*params.dydx(1);

% Optional output
varns={f1,f2,x,y};
varargout=varns(1:nargout);

if nargout==0
  % Plot if no output requested
  clf
  imagesc(x,y,f1)
  colormap bone
  axis image; xlabel('x'); ylabel('y'); longticks(gca); shrink(gca,1.1,1.1)
end

