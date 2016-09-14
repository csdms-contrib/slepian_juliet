function Sk=maternos(k,th,varargin)
% Sk=MATERNOS(k,[s2 nu rho],d)
%
% Calculates the three-parameter isotropic d-dimensional Matern
% spectral density used by Olhede & Simons (2013).
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter, aka sigma^2 
%          th(2)=nu   The second Matern parameter 
%          th(3)=rho  The third Matern parameter 
% d        Dimension [default is 2]
%
% OUTPUT:
%
% Sk       A column vector with all the wavenumbers unwrapped
%
% SEE ALSO:
%
% MATERNPRC, MATERNOS2D
%
% EXAMPLE:
%
% [~,th0]=simulosl;
% difer(maternos(0,th0)-th0(1)*pi*th0(3)^2/4,8,[],NaN)
%
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2015

% These are always the last three elements of the input 
s2=th(end-2);
nu=th(end-1);
rho=th(end);

% Change dimension if you like
if nargin==3
  d=varargin{1};
else
  d=2;
end

% Adjust for dimensionality by specialization
switch d
 case 2
  pd=nu;
 otherwise
  pd=gamma(nu+d/2)/gamma(nu);
end

% Calculate the denominator in the spectral density
avark=4*nu/pi^2/rho^2+k(:).^2;
% Calculate the d-dimensional spectral density
Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rho)^(2*nu).*avark.^(-nu-d/2);

