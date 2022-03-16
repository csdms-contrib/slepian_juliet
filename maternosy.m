function Cy=maternosy(y,th,varargin)
% Cy=MATERNOSY(y,th,d)
%
% Calculates the three-parameter isotropic d-dimensional Matern covariance
% used by Olhede & Simons (2013). Independent of d. Fully vectorized.
%
% INPUT:
%
% y        Lag parameter, the distance between spatial positions,
%          e.g. from XXPDIST or SSPDIST [m], can be any dimension
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2    The first Matern parameter [variance in units^2]
%          nu    The second Matern parameter [differentiability]
%          rh    The third Matern parameter [range in m]
% d        The dimensionality. Not used, for symmetry with MATERNOS only 
%
% OUTPUT:
%
% Cy       The spatial Matern covariance at all the requested lags
%
% SEE ALSO:
%
% MATERNOS, MATERNOSP, MATERNPRC
%
% EXAMPLE:
%
% [Hx,th0,p]=simulosl; bigN=1000;
% y=linspace(0,sqrt(prod(p.dydx))*sqrt(prod(p.NyNx)),bigN);
%% Discrete approximation of the integral
%% Compared to what we think it should be:
% [sum(y.*maternosy(y,th0))*(y(2)-y(1))/(2*pi) maternos(0,th0)]
%
% Last modified by fjsimons-at-alum.mit.edu, 11/03/2016
%
%% The Fourier relation is impossible to verify... since we never observe
%% all lags fine enough... and that is the point of estimating it
%% differently. See the comparison in SIMULOSL1('demo1') which is more fruitful.
% th=[2000 0.5 1];
% p.NyNx=[4300 5500]; lY=13; lX=24;
% [k,kx,ky,dci]=knum2(p.NyNx,[lY lX]); 
% p.dydx=[lY/(p.NyNx(1)-1) lX/(p.NyNx(2)-1)];
% x=[-floor(p.NyNx(2)/2):1:+floor(p.NyNx(2)/2)-1]*p.dydx(2);
% y=[-floor(p.NyNx(1)/2):1:+floor(p.NyNx(1)/2)-1]*p.dydx(1);
% [X,Y]=meshgrid(x,y); yy=sqrt(X.^2+Y.^2);
% Sb=maternos(k,th,2); Sbb=v2s(Sb,p);
% Cy=maternosy(yy,th,2); difer(Cy(dci(1),dci(2))-th(1))
% Sk=tospace(Cy,p); Skk=fftshift(v2s(Sk,p));
% plot(log10(Sbb(dci(1),:))); hold on ; plot(log10(abs(Skk(dci(1),:))));
% hold off

% These are always the last three elements of the input 
s2=th(end-2);
nu=th(end-1);
rh=th(end  );

% The argument, make sure it is a distance
argu=2*sqrt(nu)/pi/rh*abs(y);
% The evaluation
Cy=2^(1-nu)*s2/gamma(nu)*argu.^nu.*besselk(nu,argu);
% Supply the smallest arguments
Cy(y==0)=s2;
