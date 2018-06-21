function Hk=tospec(Hx,p)
% Hk=TOSPEC(Hx,p)
%
% Converts a vector of spatial-domain quantities to the physical spectral
% domain, with unitary normalization, and units, returned as a vector.
%
% INPUT:
%
% Hx    Real-valued column vector of unwrapped spatial-domain quantities 
% p     A parameter structure, containing, at a minimum,
%       NyNx  The number of samples in the y and x directions 
%       dydx  The sampling interval in the y and x directions [m m]
%
% OUTPUT:
%
% Hk     Complex-valued column vector of unwrapped Fourier-domain entries  
%
% SEE ALSO:
%
% TOSPACE, KNUMS, KNUM2
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/20/2018

% Unitary transform, and it has physical units
Hk=indeks(fftshift(fft2(reshape(Hx,p.NyNx))),':')...
   /sqrt(prod(p.NyNx))*sqrt(prod(p.dydx));
