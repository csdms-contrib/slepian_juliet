function Cy=maternosy(y,th,varargin)
% Cy=MATERNOSY(y,th,d)
%
% Calculates the three-parameter isotropic d-dimensional Matern
% covariance used by Olhede & Simons (2013). Independent of d. 
%
% INPUT:
%
% y        Lag parameter, the distance between spatial positions,
%          e.g. from XXPDIST or SSPDIST [m]
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
% Last modified by fjsimons-at-alum.mit.edu, 10/31/2016

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

