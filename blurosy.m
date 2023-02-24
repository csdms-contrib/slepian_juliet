function [Sbar,k]=blurosy(th,params,xver,method)
% [Sbar,k]=blurosy(th,params,xver,method)
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
%         'efs' exact, efficient and faster, exploiting Hermiticity [default]
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
defval('xver',1)
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
    [Cyy,t]=spatmat(ycol,xrow,th,params,xver);

    % Here is the blurred covariance on the 'double' grid
    Hh=fftshift(realize(fft2(ifftshift(Cyy))));

    % Play with a culled DFTMTX? Rather now subsample to the 'complete' grid
    Hh=Hh(1+mod(NyNx(1),2):2:end,1+mod(NyNx(2),2):2:end);
  case 'efs'
    % Fully exact and trying to be faster for advanced symmetry in the covariance
    ycol=[0:NyNx(1)-1]';
    xrow=[0:NyNx(2)-1] ;
    
    % Here is the Matern spatial covariance on the distance grid again, see above
    [Cyy,t]=spatmat(ycol,xrow,th,params,xver);

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
function [Cyy,t]=spatmat(ycol,xrow,th,params,xver)
% Returns the Modified spatial covariance whose Fourier transform is
% the blurred spectrum after spatial data tapering, i.e. the expected
% periodogram. Also returns the autocovariance of the applied taper in
% case it was not a boxcar.

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
    keyboard
    tx=params.taperx;
    txx=zeros(size(tx));
    % Produce the normalized autocorrelation sequence eq. (12)
    for i=1:NyNx(1)
        for j=1:NyNx(2)
            txx(i,j)=sum(sum((tx(1:NyNx(1)-i+1,1:NyNx(2)-j+1)).*(conj(tx(i:end,j:end)))));
        end
    end
    % Normalize - that probably should be using the values and not just
    % counting them, as is of course the case for a unit window
    txx=txx/prod(NyNx);

    % Here too should use FFT where we can, see compute_kernels
    % internally and below
    if xver==1
        
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

    if xver==1
        % Try alternative for the boxcar all-ones
        tx=ones(NyNx)/sqrt(prod(NyNx));
        % Need to cut one off Arthur says, possibly need to
        % re-re-visit these even/odd comparisons in BLUROS, if it ever
        % gets to that point; currently the comparison is favorable
        t2=fftshift(ifft2(abs(fft2(tx,2*size(tx,1)-1,2*size(tx,2)-1)).^2));
        % Fix the rim by adding zeroes top and left
        t2=[zeros(size(t2,1)+1,1) [zeros(1,size(t2,2)) ; t2]];
        % Check the difference between these two implementations,
        % all checked for even/odd/method combinations on 2/24/2023
        if size(t)==NyNx
            % Then the doubling inside this block needs to be undone
            t2=t2(NyNx(1)+1:end,NyNx(2)+1:end);
        end
        diferm(t,t2); disp('I checked')
    end
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
