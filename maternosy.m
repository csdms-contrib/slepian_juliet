function Cy=maternosy(y,th0s,varargin)
% Cy=MATERNOSY(y,th0s,d)
%
% Calculates the three-parameter isotropic d-dimensional Matern
% covariance used by Olhede & Simons (2013). Independent of d. 
%
% INPUT:
%
% y        Lag parameter, the distance between spatial positions
% th0s     The spectral parameter vector, containing
%          s2    The first Matern parameter, aka sigma^2
%          nu    The second Matern parameter
%          rh    The third Matern parameter
% d        The dimensionality. Leave it out, for symmetry only
%
% OUTPUT:
%
% Cy       A column vector with all the wavenumbers unwrapped
%
% SEE ALSO:
%
% MATERNOS, MATERNOS2D
%
% EXAMPLE:
%
% [Hx,th0,p]=simulosl; bigN=1000;
% y=linspace(0,sqrt(prod(p.dydx))*sqrt(prod(p.NyNx)),bigN);
%% Discrete approximation of the integral
%% Compared to what we think it should be:
% [sum(y.*maternosy(y,th0))*(y(2)-y(1))/(2*pi) maternos(0,th0)]
%
% Last modified by fjsimons-at-alum.mit.edu, 07/15/2015

% These are always the last three elements of the input 
s2=th0s(end-2);
nu=th0s(end-1);
rho=th0s(end);

% t=tic;
% The argument, make sure it is a distance
argu=2*sqrt(nu)/pi/rho*abs(y);
% The evaluation
Cy=2^(1-nu)*s2/gamma(nu)*argu.^nu.*besselk(nu,argu);
% Supply the smallest arguments
Cy(y==0)=s2;

% disp(sprintf('%s took %f seconds',upper(mfilename),toc(t)))
