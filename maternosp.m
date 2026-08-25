function [Sk,k]=maternosp(th,params,xver,tsto)
% [Sk,k]=maternosp(th,params,xver,tsto)
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
%               -1 Space multiplicative BLUROSY using the exact procedure
%              Inf Rather simulate using SGP invariant embedding, error
%           taper  0 there is no taper near of far
%                  1 it's a unit taper, implicitly
%                  OR an appropriately sized taper with proper values 
%                    (1 is yes and 0 is no and everything in between)
% xver     Perform excessive verification [0, 1 or 2]
% tsto     Offset in Fourier indices to be passed on to BLUROSY
%
% OUTPUT:
%
% Sk       A column vector with all the wavenumbers unwrapped, OR
%          a matrix if it is the multivariable lithospheric case
% k        The wavenumber matrix (the norm of the wave vectors), unwrapped
% 
% SEE ALSO:
%
% MATERNOS, MATERNOSY, MATERNPRC, KNUMS, BLUROS, BLUROSY, BLURCHECK
%
% Last modified by owalbert-at-princeton.edu, 03/18/2025
% Last modified by fjsimons-at-alum.mit.edu, 06/25/2026

% Default is to oververify it all
defval('xver','1')

% We only need to specify tsto for some variance calculations methods
defval('tsto',[]);

% The dimensionality of the problem (not being actively used)
d=2;

% Extract the needed parameters of the simulation variables
blurs=params.blurs;
try
    % Make sure the lithospheric case is testable
    DEL=params.DEL;
end

%  Calculate the Matern spectrum for the spectral parameters
switch blurs
  case {0,1}
    % The unblurred-property wavenumber grid
    k=knums(params); k=k(:);
    % disp(sprintf('%s without blurring',upper(mfilename)))
    Sk=maternos(k,th,[],d);
    % If it's the lithospheric case
    if exist('DEL','var')
        % The Cholesky decomposition of the lithospheric-spectral matrix
        [~,~,~,T]=Tros(k,th,params); 
        % The spectral matrix in case you might want it
        Sk=[Sk.*T(:,1) Sk.*T(:,2) Sk.*T(:,3)];
    end
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
       
        % If it's the lithospheric case
        if exist('DEL','var')
            % Something else is going to need to happen to make a "new way"
            % incorporate the "oldway" of SIMULROS
            
            % The lithospheric-spectral matrix on this second grid
            [~,~,~,T]=Tros(knums(params,1),th,params); 

            % FJS should make the gravity BEFORE blurring says SCO
            % SG2=[2*pi*G*DEL(2)*exp(-k(:).*z2)]   .*S(:,2);
            % SG3=[2*pi*G*DEL(2)*exp(-k(:).*z2)].^2.*S(:,3);

            % Note this gets explicitly expanded
            [Sk,k]=bluros(maternos(knums(params,1),th,[],d).*T,params,xver);
        else
            % Only for the BOXCAR at this point
            [Sk,k]=bluros(maternos(knums(params,1),th,[],d),params,xver);
        end
    else
        % Here is the alternative EXACT way of doing it, which does away
        % with the approximate convolutional refinement procedure
        % disp(sprintf('%s with exact blurring',upper(mfilename)))
        [Sk,k]=blurosy(th,params,xver,[],tsto);
    end
end
