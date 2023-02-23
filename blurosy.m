function [Sbar,k]=blurosy(th,params,xver,method)
% [Sbar,k]=blurosy(th,params,xver,method)
%
% Wavenumber blurring of a univariate Matern covariance structure with
% the periodogram of a boxcar. This is the exact (fast) explicit way
% which requires no convolutional grid refinement. Later, will build
% in other types of windows. Note this is the expected
% periodogram. Unlike BLUROS, this IS a stand-alone code, in that it
% is not blurring an input parameterized-spectral-density but rather
% producing a blurred-spectral-density from input spectral-parameters.
%
% INPUT:
%
% th      The spectral parameter vector with elements:
%         th(1)=s2   The first Matern parameter, aka sigma^2 
%         th(2)=nu   The second Matern parameter 
%         th(3)=rho  The third Matern parameter 
% params  Parameters of this experiment, the ones that are needed are:
%         dydx  sampling interval in the y and x directions [m m]
%         NyNx  number of samples in the y and x directions
%         blurs -1 as appropriate for this procedure, any other errors
%         taperx  0 there is no taper near of far
%                 1 it's a unit taper, implicitly
%                OR an appropriately sized taper with proper values
% xver    1 Extra verification via BLURCHECK
%         0 No checking at all
% method  'ef' exact, efficient and fast
%         'efs' exact, efficient and faster, exploiting Hermiticity [default]
%         'tocome' Arthurs' latest, after January 2023 discussion
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, on the original requested
%         dimension as identified by 'params' from the input
% k       The wavenumber matrix (the norm of the wave vectors), unwrapped
%
% SEE ALSO:
%
% SIMULOSL, BLUROS, MATERNOSP, BLURCHECK
%
% EXAMPLE:
%
% BLUROS('demo1',pp,bb) compares against BLUROS where
%                pp     A square matrix size (one number)
%                   bb  A MATERNOSP blurring densification (one number)
%
% Last modified by arthur.guillaumin.14-at-ucl.ac.uk, 10/15/2017
% Last modified by fjsimons-at-alum.mit.edu, 02/23/2023

if params.blurs>=0 & ~isinf(params.blurs)
  error('Are you sure you should be running BLUROSY, not BLUROS?')
end
if isinf(params.blurs)
  error('Are you sure you should be running BLUROSY, not MATERNOSY?')
end

% Defaults
defval('xver','1')
% It is here and now only that we decide to always go with 'efs'
defval('method','efs')

% Target dimensions, the original ones
NyNx=params.NyNx;
dydx=params.dydx;

switch method 
 case 'ef'
  % http://blogs.mathworks.com/steve/2010/07/16/complex-surprises-from-fft/

  % Fully exact and not particularly fast, still much faster than BLUROS
  ycol=[-NyNx(1):NyNx(1)-1]';
  xrow=[-NyNx(2):NyNx(2)-1] ;

  % Here is the Matern spatial covariance on the distance grid,
  % multiplied by the spatial taper in a way that its Fourier
  % transform can be the convolution of the spectral density with the
  % spectral density of the taper, i.e. the expected periodogram
  Cyy=spatmat(ycol,xrow,th,params);

  % Here is the blurred covariance on the 'double' grid
  Hh=fftshift(realize(fft2(ifftshift(Cyy))));

  % Play with a culled DFTMTX? Rather now subsample to the 'complete' grid
  Hh=Hh(1+mod(NyNx(1),2):2:end,1+mod(NyNx(2),2):2:end);
  case 'efs'
  % Fully exact and trying to be faster for advanced symmetry in the covariance
  ycol=[0:NyNx(1)-1]';
  xrow=[0:NyNx(2)-1] ;

  % Here is the Matern spatial covariance on the distance grid again, see above
  Cyy=spatmat(ycol,xrow,th,params);

  % Exploit the symmetry just a tad, which allows us to work with smaller matrices
  q1=fft2(Cyy);
  q4=q1+[q1(:,1) fliplr(q1(:,2:end))];

  % Here is the blurred covariance on the 'complete' grid
  Hh=fftshift(2*real(q4-repmat(fft(Cyy(:,1)),1,NyNx(2)))...
	      -repmat(2*real(fft(Cyy(1,1:end))),NyNx(1),1)...
	      +Cyy(1,1));
end

% Normalize and vectorize
Sbar=Hh(:)*prod(dydx)/(2*pi)^2;

% Check Hermiticity of the result
if xver==1
  blurcheck(Sbar,params)
end

% Produce the unwrapped wavenumbers if you've requested them to be output
if nargout>1
  k=knums(params);
  k=k(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cyy,t]=spatmat(ycol,xrow,th,params)
% Returns the modified spatial covariance whose transform is the
% blurred spectrum when the taper is a rectangular boxcar, i.e. the
% expected periodogram, and also returns the applied taper

% Target dimensions, the original ones
NyNx=params.NyNx;
dydx=params.dydx;

% Here is the distance grid
y=sqrt(bsxfun(@plus,[ycol*dydx(1)].^2,[xrow*dydx(2)].^2));

% Specify the spatial taper
if prod(size(params.taperx))>1
    % Completely general windows
    % See Arthur's note for more general windows, use IFF2/FFT2 you need, see
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/NonParametricEstimation/Periodogram.m
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/QuickRunX/expected_periodogram.m
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/QuickRunX/kernel_modulation.m

else
    % It's a unit spatial taper operation, the triangles coming out of
    % the convolution of the unit window functions, the triangle c_g,n(u) of
    % eqs 12-13 in Guillaumin et al, 2022, doi: 10.1111/rssb.12539
    triy=1-abs(ycol)/NyNx(1);
    trix=1-abs(xrow)/NyNx(2);
    % Here is the triangle grid c_g,n(u) of 
    t=bsxfun(@times,triy,trix);

    % % Try alternative for the boxcar all-ones
    % I=ones(NyNx)/sqrt(prod(NyNx));
    % % Need to cut one off Arthur says
    % t2=fftshift(ifft2(abs(fft2(I,2*size(I,1)-1,2*size(I,2)-1)).^2));
    % % Fix the rim by adding zeroes top and left
    % t2=[zeros(size(t2,1)+1,1) [zeros(1,size(t2,2)) ; t2]];

end

% Need the modified spatial covariance
Cyy=maternosy(y,th).*t;

% Arthur: diag(U* * Cyy * U) where U is the DFMTX
% is the expected periodogram but to get the variance of the gradient we need the whole thing
% (U* * Cyy * U)
% of course in two dimensions
% periodogram is diag of covariance of dft

% keyboard

% % I would be weighted in case of irregular sampling

% % Generic mask using blobi
% %I=blobi;


% diferm(t,t2)
% clf

% subplot(221)
% imagesc(t); axis image; title('t')
% subplot(222)
% imagesc(t2); axis image; title('t2')
% % No sense in comparing the rims
% subplot(223)
% imagesc(t(2:end,2:end)-t2(2:end,2:end)); axis image; title('difference'); colorbar
% subplot(224)
% imagesc(t(2:end,2:end)./t2(2:end,2:end)); axis image; title('ratio'); colorbar
% caxis(1+[-1 1]*1e-10)

% % A third way for a quadrant
% for nrow=1:NyNx(1)
%     for ncol=1:NyNx(2)
%         % Autocorrelation of the index sequence, i.e. spatial version of above
%         t3(nrow,ncol)=sum(sum((I(1:NyNx(1)-nrow+1,1:NyNx(2)-ncol+1)).*(conj(I(nrow:end, ncol:end)))));
%     end
% end
% % t3=t3/sqrt(prod(NyNx));
% diferm(t(51:100,51:100),t3)
