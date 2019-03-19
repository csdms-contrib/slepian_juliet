function [k,kx,dci,dcn,dx]=knum1(n,pl)
% [k,kx,dci,dcn,dx]=knum1(n,lX)
%
% Angular wavenumber (not frequency!) axis for a 1D Fourier transform,
% suitable after FFTSHIFT and for manipulation ensuring real-field Hermitian
% symmetry under IFFT. The wavenumbers increase symmetrically away from
% zero, but the Nyquist frequency is only reached at the left for an even
% number of samples.
%
% INPUT:
% 
% n         Sample length of FFT vector, equal to length of original vector 
% lX        Physical lenggth of original matrix
%
% OUTPUT:
%
% k         The wavenumber vector (the absolute values of the wave vectors)
% kx        The components of the wavenumber vector
% dci       The index to the DC component [at floor(n/2)+1]
% dcn       The index to the components at exactly the Nyquist
% dx        The sampling interval in the space domain
%
% The "corner Nyquist" rate is at 1 and the DC rate is the "center
% point", dci, which is on the left quadrant for even lengths; the way
% IFFT/FFTSHIFT expect it. The negative x frequencies are on the left.
%
% SEE ALSO:
%
% KNUM2, the equivalent two-dimensional function
%
% EXAMPLE:
%
% N=60+round(rand); h=rindeks(peaks(N),randi(N)); H=fftshift(fft(h)); F=randn*10;
% hf=ifft(ifftshift(H.*exp(-F*knum1(N,100)))); plot(hf)
% 
% EXAMPLE:
%
% H=rindeks(peaks(200),randi(200)); Hk=fftshift(fft(H));
% [K,kx,dci,dcn,dx]=knum1(length(Hk),length(H)-1);
%% Where and what is the "DC" component?
% difer(sum(H)-Hk(dci))
%% How is Parseval's theorem satisfied?
% difer(H(:)'*H(:)-Hk(:)'*Hk(:)/length(Hk),8)
%
%% The equivalence of FFTAXIS1D and KNUM1 is explicit:
% N=256;
% [fax,selekt]=fftaxis1D(rand(N,1),N,N-1);
% [k,kx]=knum1(N,N-1);
% fx=-fliplr(indeks(kx,selekt)/2/pi);
% difer(fx(:)-fax(:))
%
%% This is only subtly different from FFTAXIS, which is becoming obsolete:
% [xfaks,fnx,xsint]=fftaxis(length(H),length(Hk),length(H)-1);
% difer(kx/2/pi+fliplr(xfaks))
%
% SEE ALSO: FFTAXIS, FFTAXIS1D, KNUM2, KNUMS, RANDGPN, BRACEWELL
%
% Last modified by fjsimons-at-alum.mit.edu, 03/23/2018

% Supply defaults for testing
defval('n',32+round(rand))
defval('pl',n-1)

% Vector size
N=n;

% Physical dimension
lx=pl;

% Sampling interval
dx=lx/(N-1);

% Wave vector axis
kx=2*pi*linspace(-floor(N/2),floor((N-1)/2),N)/N/dx;

% Wavenumbers
k=abs(kx);

% The index of the zero frequency, guaranteed
dci=floor(N/2)+1;
difer(k(dci),[],[],NaN)

% The index where the Nyquist is exactly reached, which is only the
% left side if the dimension is even.
dcn=1;
isev=~mod(n,2);
if isev
  difer(k(dcn)/2/pi*pl/(n-1)-1/2,[],[],NaN)
end
dcn=dcn(isev);
