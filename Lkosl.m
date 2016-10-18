function [Lk,Xk]=Lkosl(k,th,params,Hk)
% [Lk,Xk]=Lkosl(k,th,params,Hk)
%
% Computes the likelihood function for the SINGLE-FIELD isotropic Matern
% model in Olhede & Simons (2013). Takes out the zero wavenumber. 
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
%          The zero wavenumber is present but excluded from consideration
% th       The parameter vector with elements [unscaled]:
%          s2  The first Matern parameter, aka sigma^2
%          nu  The second Matern parameter
%          rho The third Matern parameter
% params   A structure with AT LEAST these constants that are known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
% Hk       A [prod(params.NyNx)*1]-column vector of complex Fourier-domain observations
%
% OUTPUT:
%
% Lk       A one-column likelihood vector with the wavenumbers unwrapped
% Xk       A quadratic piece of it that gets used in the analysis of residuals
%
% SEE ALSO: 
%
% LOGLIOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 10/18/2016

% Extra verification?
defval('xver',1)

% Extract the needed parameters of the simulation variables
blurs=params.blurs;
NyNx=params.NyNx;

% We need the power spectrum and its ratio to the observations 

switch blurs
 case {0,1}
  % Too much screen time, FMINUNC/FMINCON evaluate this a lot
  % disp(sprintf('%s without blurring',upper(mfilename)))
  %  Calculate the Matern spectrum with the spectral parameters
  S=maternos(k,th);
 otherwise
  % Too much screen time, FMINUNC/FMINCON evaluate this a lot
  % disp(sprintf('%s with blurring factor %i',upper(mfilename),blurs))
  if blurs>1
    % Blurs IS the refinement parameter; make new wavenumber grid
    % Now make the spectral-spectral portion of the spectral matrix
    % Which we need to convolve now in two dimensions
    % And then do subsampling onto the original target grid
    S=bluros(maternos(knums(params,1),th),params);
  else
    % Here is the alternative EXACT way of doing it, which does away
    % with the approximate convolutional refinement procedure
    S=blurosy(th,params);
  end
end

% The average of Xk needs to be close to one as will be tested 
warning off MATLAB:log:logOfZero
Xk=abs(Hk).^2./S;
warning on MATLAB:log:logOfZero
Lk=-log(Sb)-Xk;

% Trouble is at the central wave numbers, we take those out 
% Find the zero wavenumber
kzero=sub2ind(NyNx,floor(NyNx(1)/2)+1,floor(NyNx(2)/2)+1);
% Check that this is correctly done
difer(k(kzero),[],[],NaN)
% Behavior is rather different if this is NOT done... knowing that it
% will not blow up but rather be some numerically large value
Lk(kzero)=NaN;

% Should make sure that this is real! Don't take any chances
Lk=realize(Lk);
