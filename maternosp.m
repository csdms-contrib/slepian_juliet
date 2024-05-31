function [Sk,k]=maternosp(th,params,xver)
% [Sk,k]=maternosp(th,params,xver)
%
% Calculates the three-parameter isotropic d-dimensional Matern spectral
% density used by Olhede & Simons (2013), with or without blurring, exact or
% not. Wavenumber grid is complete and depends on your input parameters.
%
% INPUT:
%
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2   The first Matern parameter [variance in units^2]
%          nu   The second Matern parameter [differentiability]
%          rho  The third Matern parameter [range in m]
% params   A structure with AT LEAST these constants that are known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 No wavenumber blurring
%                1 No wavenumber blurring, effectively
%                N Fejer convolutional  BLUROS  on an N-times refined grid
%               -1 Fejer multiplicative BLUROSY using the exact procedure
%              Inf Rather simulate using SGP invariant embedding, error
%           taper  0 there is no taper near of far
%                  1 it's a unit taper, implicitly
%                  OR an appropriately sized taper with proper values 
%                    (1 is yes and 0 is no and everything in between)
% xver     Perform excessive verification [0, 1 or 2]
%
% OUTPUT:
%
% Sk       A column vector with all the wavenumbers unwrapped
% k        The wavenumber matrix (the norm of the wave vectors), unwrapped
% 
% SEE ALSO:
%
% MATERNOS, MATERNOSY, MATERNPRC, KNUMS, BLUROS, BLUROSY, BLURCHECK
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2023

% Default is to oververify it all
defval('xver','1')

% The dimensionality of the problem (not being actively used)
d=2;

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

%  Calculate the Matern spectrum for the spectral parameters
switch blurs
  case {0,1}
    % The unblurred-property wavenumber grid
    k=knums(params); k=k(:);
    % disp(sprintf('%s without blurring',upper(mfilename)))
    Sk=maternos(k,th,d);
  otherwise
    if isinf(blurs)
        error('You should be calling MATERNOSY, not MATERNOSP!')
    end
    if blurs>1 & ~isinf(blurs)
        % disp(sprintf('%s with blurring factor %i',upper(mfilename),blurs))

        % If I stay on the same k-grid, I'm really not doing any convolution at
        % all as you can see quickly. So by "suitably discretizing" the
        % convolutional operator we mean: performing it on a highly densified
        % grid of which the target grid may be a subset. Doing this on the
        % same grid would be the "inverse crime" of NOT changing the grid
        % all. Run Fk for this case to see it then would be a delta function.
        
        % Blurs IS the refinement parameter; make new wavenumber grid. Make the
        % spectral-spectral portion of the spectral matrix, convolve in two
        % dimensions, and then subsample onto the original target grid

        % Only for the BOXCAR at this point
        [Sk,k]=bluros(maternos(knums(params,1),th),params,xver);
    else
        % Here is the alternative EXACT way of doing it, which does away
        % with the approximate convolutional refinement procedure
        % disp(sprintf('%s with exact blurring',upper(mfilename)))
        [Sk,k]=blurosy(th,params,xver);
    end
end
