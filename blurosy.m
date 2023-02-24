function [Sbar,k,t]=blurosy(th,params,xver,method)
% [Sbar,k,t]=blurosy(th,params,xver,method)
%
% Wavenumber blurring of a univariate Matern spectral density with the
% periodogram of a spatial taper. This is the exact (fast) explicit
% way which requires no convolutional grid refinement. Note that the
% result is the expected periodogram. Unlike BLUROS, this IS a
% stand-alone code, in that it is not blurring an input parameterized
% spectral density but rather producing a blurred spectral density
% from input spectral-parameters, through Fourier transformation of
% the product of the Matern correlation function with the
% autocorrelation of the taper window function. Equations refer to
% Guillaumin et al., 2022, doi: 10.1111/rssb.12539
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
%         'efs' exact, efficient and faster, exploiting symmetry [default]
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, on the original requested
%         dimension as identified by 'params' from the input
% k       The wavenumber matrix (the norm of the wave vectors), unwrapped
% t       The autocorrelation of the spatial taper... 
%         may never need it explicitly, used in SIMULOSL and LOGLIOSL
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
defval('xver',1)
defval('method','efs')

% Target dimensions, the original ones
NyNx=params.NyNx;
dydx=params.dydx;

switch method 
  case 'ef'
    % Generates a double grid from which we subsample
    % Fully exact and not particularly fast, still much faster than BLUROS

    % Here are the full lags
    ycol=[-NyNx(1):NyNx(1)-1]';
    xrow=[-NyNx(2):NyNx(2)-1] ;

    % Here is the Matern spatial covariance on the double distance grid,
    % multiplied by the spatial taper in a way that its Fourier
    % transform can be the convolution of the spectral density with the
    % spectral density of the taper, i.e. the expected periodogram
    [Cyy,tyy]=spatmat(ycol,xrow,th,params,xver);

    % http://blogs.mathworks.com/steve/2010/07/16/complex-surprises-from-fft/
    % Here is the blurred covariance on the 'double' grid
    Hh=fftshift(realize(fft2(ifftshift(Cyy))));

    % Play with a culled DFTMTX? Rather now subsample to the 'complete' grid
    Hh=Hh(1+mod(NyNx(1),2):2:end,1+mod(NyNx(2),2):2:end);
  case 'efs'
    % Generates a sample-size grid by working from a quarter, rest symmetric
    % Fully exact and trying to be faster for advanced symmetry in the covariance

    % Here are the partial lags
    ycol=[0:NyNx(1)-1]';
    xrow=[0:NyNx(2)-1] ;
    
    % Here is the Matern spatial covariance on the quarter distance grid, see above
    [Cyy,tyy]=spatmat(ycol,xrow,th,params,xver);

    % Exploit the symmetry just a tad, which allows us to work with smaller matrices
    q1=fft2(Cyy);
    q4=q1+[q1(:,1) fliplr(q1(:,2:end))];

    % Here is the blurred covariance on the 'complete' grid
    Hh=fftshift(2*real(q4-repmat(fft(Cyy(:,1)),1,NyNx(2)))...
	        -repmat(2*real(fft(Cyy(1,1:end))),NyNx(1),1)...
	        +Cyy(1,1));
    % If you ever wanted t to come out you'll need to unquarter it
    % t=
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
function [Cyy,t]=spatmat(ycol,xrow,th,params,xver)
% [Cyy,t]=spatmat(ycol,xrow,th,params,xver)
%
% Returns the modified spatial covariance whose Fourier transform is
% the blurred spectrum after spatial data tapering, i.e. the expected
% periodogram. Also returns the autocovariance of the applied taper
% which is an essential part of this operation

% Dimensions of the original grid
NyNx=params.NyNx;
dydx=params.dydx;

disp('Next force us in the loop even for unit taper')

% Specify the spatial taper
kb
if prod(size(params.taperx))>1
    % Now compare with the other mechanisms
    % Completely general windows
    % See Arthur's note for more general windows, use IFF2/FFT2 you need, see
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/NonParametricEstimation/Periodogram.m
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/QuickRunX/expected_periodogram.m
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/QuickRunX/kernel_modulation.m
    tx=params.taperx;
    % Make sure it's normalized
    tx=tx./sum(sum(tx.^2));

    % Produce the normalized autocorrelation sequence eq. (12)
    t=zeros(size(tx));
    % This is directly compatible with the efs but for ef will need
    % to do something else, however, no rocket science
    %    for i=1:NyNx(1)
    %       for j=1:NyNx(2)
    % It's quite vital that these be colon ranges (faster) or (like
    % here) row index vectors... row/column won't work
    for i=xrow(:)'+1
        for j=ycol(:)'+1
            % Vectorize? Check out XCORR2, that's good
            t(i,j)=sum(sum(tx(1:NyNx(1)-i+1,1:NyNx(2)-j+1)).*(conj(tx(i:end,j:end)))));
        end
    end
    % Normalize - that probably should be using the values and not just
    % counting them, as is of course the case for a unit window
    t=t/prod(NyNx);

    % Here too should use FFT where we can, see compute_kernels
    % internally and below, I would image that's just the same thing
    if xver==1
        disp('Do something or skip')
    end
else
    % Just a single number, could be 0 or 1 both telling us "not" spatially
    % tapered which is of course the same as a "unit" taper, in essence
    % It's a unit spatial taper operation, the triangles coming out of
    % the autocorrelation of the unit window functions, the triangle c_g,n(u) of
    % eqs 12-13 in Guillaumin et al, 2022, doi: 10.1111/rssb.12539
    triy=1-abs(ycol)/NyNx(1);
    trix=1-abs(xrow)/NyNx(2);
    % Here is the gridded triangle for this case
    t=bsxfun(@times,triy,trix);

    % It was a unit taper all along
    if xver==1
        tx=ones(NyNx)/sqrt(prod(NyNx));
    end
end

if xver==1
    % Need to cut one off Arthur says, possibly need to
    % re-re-visit these even/odd comparisons in BLUROS, if it ever
    % gets to that point; currently the comparison is favorable
    t2=fftshift(ifft2(abs(fft2(tx,2*size(tx,1)-1,2*size(tx,2)-1)).^2));
    % Fix the rim by adding zeroes top and left
    t2=[zeros(size(t2,1)+1,1) [zeros(1,size(t2,2)) ; t2]];
    % Check the difference between these two implementations,
    % all checked for even/odd/method combinations on 2/24/2023
    if all(size(t)==NyNx)
        % Then the doubling inside this block needs to be undone
        t2=t2(NyNx(1)+1:end,NyNx(2)+1:end);
    end
    diferm(t,t2); %disp('I checked')
end

% Here is the distance grid, whose size depends on the input
y=sqrt(bsxfun(@plus,[ycol*dydx(1)].^2,[xrow*dydx(2)].^2));

% Modified the spatial covariance
Cyy=maternosy(y,th).*t;

% Remind me:
% Arthur: diag(U* * Cyy * U) where U is the DFMTX
% is the expected periodogram but to get the variance of the gradient we need the whole thing
% (U * Cyy * U)
% of course in two dimensions
% periodogram is diag of covariance of dft

% keyboard

