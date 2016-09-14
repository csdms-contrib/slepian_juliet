function Sbar=bluros(S,params,xver)
% Sbar=BLUROS(S,params,xver)
%
% Spectral blurring with periodogram of a boxcar. If we're talking about a
% matrix, there are three columns, and there can be an eigenvalue check.
% This is the approximate, slower, convolutional way which requires a "grid
% refinement" parameter. Later, will build in other types of windows.
%
% INPUT:
%
% S       The spectral matrix with its wavenumbers unwrapped, on a
%         wavenumber grid that was refined from the original
%         For univariate spectra, it's a column vector. 
%         For bivariate spectra, it's [S_XX(k(:)) SXY(k(:)) SYY(k(:))]
% params  Parameters of this experiment, the ones that are needed are:
%         dydx  sampling interval in the y and x directions [m m]
%         NyNx  number of samples in the y and x directions
%         blurs 0 Don't blur likelihood using the Fejer window
%               1 With no refinement, this is like not blurring at all!
%               N Blur likelihood using the Fejer window [default: N=2]
% xver    1 Extra verification, among other ways, via BLURCHECK
%         0 No checking at all
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, interpolated to original requested
%         dimension as identified by 'params' from the input
%
% SEE ALSO:
%
% SIMULOSL, BLUROSY
%
% NOTES:
%
% It looks to me that keeping blurs to an ODD number performs the best from
% either even- or odd-sized inputs, despite the parity fix effectuated in
% KNUMS. This from running and plotting lots of simulations in EGGERS1.
%
% Maybe should formally increase it in those cases so as to never worry?
%
% Last modified by fjsimons-at-alum.mit.edu, 09/13/2016

if params.blurs<0
  error('You should be running BLUROSY, not BLUROS!')
end

% Set defaults
defval('xver',1)

% Target dimensions, the original ones
NyNx=params.NyNx;
dydx=params.dydx;

if xver==1
  % The unblurred wavenumbers
  [k,dci,~,kx,ky]=knums(params);
  % The blurred wavenumbers
  [~,kzero,~,kx2,ky2]=knums(params,1);
  % Do the check, don't bother with dcn2 as that's done inside knum2
  diferm(k(dci(1),dci(2)))
  diferm(k(kzero))
else
  % The unblurred wavenumbers
  [~,~,~,kx,ky]=knums(params);
  % The blurred wavenumbers
  [~,~,~,kx2,ky2]=knums(params,1);
end

% Find out what you've actually been given, remember S worked on KNUMS
NyNx2=[length(ky2) length(kx2)];

% Use the EFFECTIVE blurs for the Fejer kernel!
blurs=sqrt(prod(NyNx2)/prod(NyNx));

if blurs~=params.blurs && xver==1
  % Could have been innocuous, from the parity protection, or more major
  disp(sprintf('Overriding blurs from %4.2f to %4.2f from KNUMS parity',...
               params.blurs,blurs))
end

% This is the periodogram of the boxcar for the old grid 
% interpolated to the new grid so we don't just hit the nodes
% See http://blinkdagger.com/matlab/matlab-fft-and-zero-padding/
% Should do this once and save it. Now let's not forget that in our
% formalism we force the fft/ifft to be unitary
% If params.blurs were to be 1, we would get a single 1 in the center
% Note that the kernel veluas are very different depending on even/odd dimensions
Fejk=fftshift(abs(fft2(repmat(1/prod(NyNx)/blurs,NyNx),NyNx2(1),NyNx2(2))).^2);

% Make sure it is unitary and norm-preserving
difer(sum(Fejk(:))-1,[],[],NaN)

if xver==1
  % Check Hermiticity of the Fejer kernel, no reason to doubt it
  hermcheck(Fejk)
end

% Check that there is no roundoff going on in the interpolation
% This case might apply when we start from an odd number and double
if [kx(end)-kx2(end)]>0; kx2(end)=kx(end); end
if [ky(end)-ky2(end)]>0; ky2(end)=ky(end); end

% This case might apply when we start from an even number and double
if [kx(1)-kx2(1)]<0; kx2(1)=kx(1); end
if [ky(1)-ky2(1)]<0; ky2(1)=ky(1); end

% Fill to original dimensions
Sbar=nan(prod(NyNx),size(S,2));
% Make sure there are no NaNs in the output
for in=1:size(S,2)
  % Later, consider griddedInterpolant
  Hh=interp2(kx2(:)',ky2(:),...
             conv2(Fejk,reshape(S(:,in),NyNx2),'same'),...
             kx(:)',ky(:));
  Sbar(:,in)=Hh(:);
end

% Check that no extrapolation was demanded, effectively
% but know that griddedInterpolant would have EXTRApolated fine
difer(sum(isnan(Sbar(:))),[],2,NaN)

% Should use the Claerbout helix! Not so, says Sergey Fomel.
% convmtx2 needs more memory
% Actually, should look into FFTW. But also limit to halfplane.
% disp(sprintf('BLUROS %i %i %i',blurs,NyNx2(1),NyNx2(2)));

% Check Hermiticity and positive-definiteness
if xver==1
  blurcheck(Sbar,params)
end

