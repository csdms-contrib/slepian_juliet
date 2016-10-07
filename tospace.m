function Hx=tospace(Hk,p)
% Hx=TOSPACE(Hk,p)
%
% Converts a vector of spectral-domain quantities to the physical spatial
% domain, with unitary normalization, and units, returned as a vector.
%
% INPUT:
%
% Hk       Complex-valued column vector of unwrapped Fourier-domain entries  
% p        A parameter structure, containing, at a minimum
%          NyNx  The number of samples in the y and x directions
%          dydx  The sampling interval in the y and x directions [m m]
%   
%
% OUTPUT:
%
% Hx       Real-valued column vector of unwrapped spatial-domain quantities 
%
% SEE ALSO:
%
% TOSPEC, KNUMS, KNUM2
% 
% Last modified by fjsimons-at-alum.mit.edu, 10/07/2016

% Unitary transform, and it has physical units
Hx=realize(indeks(ifft2(ifftshift(reshape(Hk,p.NyNx))),':')...
	   *sqrt(prod(p.NyNx))/sqrt(prod(p.dydx)));
