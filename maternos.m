function Sk=maternos(k,th,varargin)
% Sk=MATERNOS(k,th,d)
%
% Calculates the three-parameter isotropic d-dimensional Matern
% spectral density used by Olhede & Simons (2013). Depends on d. 
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2   The first Matern parameter [variance in units^2]
%          nu   The second Matern parameter [differentiability]
%          rho  The third Matern parameter [range in m]
% d        The dimensionality [default is 2]
%
% OUTPUT:
%
% Sk       A column vector with all the wavenumbers unwrapped
%
% SEE ALSO:
%
% MATERNPRC, MATERNOSY, MATERNOSP
%
% EXAMPLE:
%
% [~,th0]=simulosl; difer(maternos(0,th0)/th0(1)-pi*th0(3)^2/4,7,[],NaN)
%
% Last modified by fjsimons-at-alum.mit.edu, 11/03/2016

% These are always the last three elements of the input 
s2=th(end-2);
nu=th(end-1);
rh=th(end  );

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
avark=4*nu/pi^2/rh^2+k(:).^2;
% Calculate the d-dimensional spectral density
Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2);

