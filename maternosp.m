function Sk=maternosp(k,th,params,xver)
% Sk=maternosp(k,th,params,xver)
%
% Calculates the three-parameter isotropic d-dimensional Matern
% spectral density used by Olhede & Simons (2013), with or without
% blurring, exact or not. 
%
%INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m], not always needed
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2   The first Matern parameter [variance in units^2]
%          nu   The second Matern parameter [differentiability]
%          rho  The third Matern parameter [range in m]
% params   A structure with AT LEAST these constants that are known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
% xver     1 Perform excessive verification
%          0 Don't 
%
% OUTPUT:
%
% Sk       A column vector with all the wavenumbers unwrapped
% 
% SEE ALSO:
%
% MATERNOS, MATERNOSY, MATERNPRC, KNUMS
%
% Last modified by fjsimons-at-alum.mit.edu, 10/31/2016

% Default is to oververify it all
defval('xver','1')

% The dimensionality of the problem (not being actively used)
d=2;

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

%  Calculate the Matern spectrum for the spectral parameters
switch blurs
 case {0,1}
  % disp(sprintf('%s without blurring',upper(mfilename)))
  Sk=maternos(k,th,d);
 otherwise
  % disp(sprintf('%s with blurring factor %i',upper(mfilename),blurs))
  if blurs>1
    % If I stay on the same k-grid, I'm really not doing any convolution at
    % all as you can see quickly. So by "suitably discretizing" the
    % convolutional operator we mean: performing it on a highly densified
    % grid of which the target grid may be a subset. Doing this on the
    % same grid would be the "inverse crime" of NOT changing the grid
    % all. Run Fk for this case to see it then would be a delta function.

    % Blurs IS the refinement parameter; make new wavenumber grid. Make the
    % spectral-spectral portion of the spectral matrix, convolve in two
    % dimensions, and then subsample onto the original target grid
    Sk=bluros(maternos(knums(params,1),th),params,xver);
  else
    % Here is the alternative EXACT way of doing it, which does away
    % with the approximate convolutional refinement procedure
    % disp(sprintf('%s with exact blurring',upper(mfilename)))
    Sk=blurosy(th,params,xver);
  end
end

