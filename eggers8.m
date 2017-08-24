function eggers8
% EGGERS7
%
% Makes FIGURES (like) 7-8 of Olhede et al. (2017), illustrating the
% chi-squaredness of the ratio residuals as did MLECHIPLOS for EGGERS6, but
% then also plots the power-spectral densities of the fields.
%
% The main engine is MLECHIPLOS
%
% Tested on 8.3.0.532 (R2014a)
%
% Last modified by fjsimons-at-alum.mit.edu, 08/24/2017

% Set parameters for creation of a data patch
fields={'dydx','NyNx','blurs','quart'};
defstruct('params',fields,{[20 20]*1e3,128*[1 1],-1,0});

% Random random parameters
th0=max(round(rand(1,3).*[1 1 4]*10),[1 1 1])./[1e-4 1 1e-4];
th0(2)=2+rand(1,1)*2;

% Create the data patch, both in spatial and Fourier domain
[Hx,~,params,k,Hk]=simulosl(th0,params); 

% Set the isotropic wavenumber cutoff to some random fraction
params.kiso=k(1,randi([params.NyNx(1)/2,params.NyNx(1)]));
% Or rather, not
params.kiso=NaN;

% Estimate the parameters via maximum-likelihood
[thhat,~,~,scl]=mleosl(Hx,[],params,[],[],[]);

% Create a title
stit=sprintf('est %s = [%i %5.2f %i]\ntru %s = [%i %5.2f %i]',...
	     '\theta',round(thhat(1)*scl(1)),thhat(2)*scl(2),round(thhat(3)*scl(3)),...
	     '\theta',round(th0(1)),th0(2),round(th0(3)));

% Produce the residuals/power spectral density figure
clf
a=mlechipsdosl(Hk,thhat,scl,params,stit);

if a==1
  % Plot the figure! EPSTOPDF doesn't do well
  disp(' ')
  figna=figdisp([],a,[],1);
  system(sprintf('ps2raster -Tf %s.eps',figna));
  system(sprintf('rm -rf %s.eps',figna));
end
