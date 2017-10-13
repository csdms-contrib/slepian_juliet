function [Sbar,k]=blurosy(th,params,xver,method)
% [Sbar,k]=blurosy(th,params,xver,method)
%
% Spectral blurring with periodogram of a boxcar, for univariate cases.
% This is the exact, fast, explicit way which requires no convolutional
% grid refinement. Later, will build in other types of windows. 
%
% INPUT:
%
% th      The true parameter vector with elements:
%         th(1)=s2   The first Matern parameter, aka sigma^2 
%         th(2)=nu   The second Matern parameter 
%         th(3)=rho  The third Matern parameter 
% params  Parameters of this experiment, the ones that are needed are:
%         dydx  sampling interval in the y and x directions [m m]
%         NyNx  number of samples in the y and x directions
% xver    1 Extra verification via BLURCHECK
%         0 No checking at all
% method  'ef' efficient and fast [default]
%         'efs' efficient and faster, due to symmetry [not working yet]
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
% Last modified by arthur.guillaumin.14-at-ucl.ac.uk, 10/15/2017
% Last modified by fjsimons-at-alum.mit.edu, 10/15/2017

if params.blurs>=0
  error('Are you sure you should be running BLUROSY, not BLUROS?')
end

% Defaults
defval('xver','1')
defval('method','ef')

% Target dimensions, the original ones
NyNx=params.NyNx;
dydx=params.dydx;

switch method 
 case 'ef'
  % http://blogs.mathworks.com/steve/2010/07/16/complex-surprises-from-fft/

  % Fully exact and not particularly fast, still faster than BLUROS
  % Distance grid, got to be careful with the incoming parity of the signal
  ycol=[-NyNx(1)+mod(NyNx(1),2):NyNx(1)-1]';
  xrow=[-NyNx(2)+mod(NyNx(2),2):NyNx(2)-1] ; 
  % The triangles coming out of the convolution... normalized?
  triy=[1-abs(ycol)/NyNx(1)];
  trix=[1-abs(xrow)/NyNx(2)];
  % Here is the distance grid
  y=sqrt(bsxfun(@plus,[ycol*dydx(1)].^2,[xrow*dydx(2)].^2));

  % Here is the triangle grid
  t=bsxfun(@times,triy,trix);
  % Here is the blurred covariance on the 'double' grid
  Hh=fftshift(realize(fft2(ifftshift(maternosy(y,th).*t))));

  % Normalize so that the C(y) is the proper transform of S(k)
  Hh=Hh*prod(dydx)/(2*pi)^2;

  % Could play with a culled DFTMTX but rather now subsample
  Hh=Hh(1:2:end,1:2:end);

  % Vectorize the argument
  Sbar=Hh(:);    
 case 'efs'
  % Fully exact and trying to be faster for advanced symmetry
  % Distance grid
  ycol=[0:NyNx(1)-1]';
  xrow=[0:NyNx(2)-1] ;
  % The triangles coming out of the convolution
  triy=1-ycol/NyNx(1);
  trix=1-xrow/NyNx(2);

  % Here is the distance grid
  y=sqrt(bsxfun(@plus,[ycol*dydx(1)].^2,[xrow*dydx(2)].^2));
  % Here is the triangle grid
  t=bsxfun(@times,triy,trix);
  
  % Need the Matern spatial covariance on the distance grid...
  % multiplied by the transform of the Fejer kernel
  sxx=maternosy(y,th).*t;

  % Pick out a column and a row
  sxc=sxx(:,1); 
  sxr=sxx(1,:);

  % Produce the exact blurred spectrum, adjust the normalization also
  [q1,q2]=deal(fft2(sxx));
  q2(:,2:end) = fliplr(q2(:,2:end));
  q4 = q1+q2-repmat(fft(sxc),1,NyNx(1));
  
  Sbar=fftshift(abs(2*real(q4)-repmat(2*real(fft(sxr))-sxx(1,1),NyNx(2),1)));
  % Normalize and vectorize
  Sbar=Sbar(:)*prod(dydx)/(2*pi)^2;
end

% Check Hermiticity
if xver==1
  blurcheck(Sbar,params)
end

% Produce the unwrapped wavenumbers if you've request them to be output
if nargout>1
  k=knums(params);
  k=k(:);
end

