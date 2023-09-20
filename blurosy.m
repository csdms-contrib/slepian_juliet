function varargout=blurosy(th,params,xver,method)
% [Sbar,k,tyy,Cyy]=blurosy(th,params,xver,method)
%
% Wavenumber blurring of a univariate Matern spectral density with the
% periodogram of a spatial taper. The result is the expected periodogram.
%
% Exact, fast, explicit way, no convolutional grid refinement (unlike BLUROS).
% Also unlike BLUROS, this IS a stand-alone code: it is not blurring an input
% parameterized spectral density but rather producing a blurred spectral density
% directly from input spectral parameters, through Fourier transformation of the
% product of the Matern correlation function with the autocorrelation of the
% taper window function.
%
% Equation numbers refer to Guillaumin et al., 2022, doi: 10.1111/rssb.12539
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
%         taper  0 there is no taper near of far
%                1 it's a unit taper, implicitly
%                OR an appropriately sized taper with explicit values 
%                   (1 is yes and 0 is no and everything in between)
% xver    1 Extra verification via BLURCHECK and alternative computations
%         0 No checking at all
% method  'ef' exact, efficient and fast
%         'efs' exact, efficient and faster, exploiting symmetry [default]
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, with the unwrapped requested dimension
%         as identified by the input 'params' (wrap with V2S)
% k       The wavenumber matrix (the norm of the wave vectors), unwrapped
% tyy     The autocorrelation of the spatial taper... which you 
%         may never need it explicitly, used in SIMULOSL and LOGLIOSL
% Cyy     The modified Matern correlation, may never need it explicitly
%
% SEE ALSO:
%
% SIMULOSL, BLUROS, MATERNOSP, BLURCHECK
%
% EXAMPLE:
%
% [H,th,p]=simulosl; p.blurs=-1;
% p.taper=1; % Implicit unit taper
% [Sbar1,k1,tyy1,Cyy1]=blurosy(th,p,1,'ef');
% [Sbar2,k2,tyy2,Cyy2]=blurosy(th,p,1,'efs');
% p.taper=ones(p.NyNx); % Explicit unit taper
% [Sbar3,k3,tyy3,Cyy3]=blurosy(th,p,1,'ef');
% [Sbar4,k4,tyy4,Cyy4]=blurosy(th,p,1,'efs');
% % All of these should be virtually identical
%
% BLUROS('demo1',pp,bb) compares against BLUROS where
%                pp     A square matrix size (one number)
%                   bb  A MATERNOSP blurring densification (one number)
%
% Last modified by arthur.guillaumin.14-at-ucl.ac.uk, 10/15/2017
% Last modified by fjsimons-at-alum.mit.edu, 09/20/2023

if params.blurs>=0 & ~isinf(params.blurs)
  error('Are you sure you should be running BLUROSY, not BLUROS?')
end
if isinf(params.blurs)
  error('Are you sure you should be running BLUROSY, not MATERNOSY?')
end

% Defaults
defval('xver',1)
defval('method','ef')

% Target dimensions, the original ones
NyNx=params.NyNx;
dydx=params.dydx;

switch method 
  case 'ef'
    % Generates a double grid from which we subsample
    % Fully exact and not particularly fast, still much faster than BLUROS

    % Here are the full lags
    ydim=[-NyNx(1):NyNx(1)-1]';
    xdim=[-NyNx(2):NyNx(2)-1] ;

    % Here is the Matern spatial covariance on the double distance grid,
    % multiplied by the spatial taper in a way that its Fourier
    % transform can be the convolution of the spectral density with the
    % spectral density of the taper, i.e. the expected periodogram
    [Cyy,tyy]=spatmat(ydim,xdim,th,params,xver);

    % http://blogs.mathworks.com/steve/2010/07/16/complex-surprises-from-fft/
    % Here is the blurred covariance on the 'double' grid
    Hh=fftshift(realize(fft2(ifftshift(Cyy))));

    % Play with a culled DFTMTX? Rather now subsample to the 'complete' grid
    Hh=Hh(1+mod(NyNx(1),2):2:end,1+mod(NyNx(2),2):2:end);
  case 'efs'
    % Generates a sample-size grid by working from a quarter, rest symmetric
    % Fully exact and trying to be faster for advanced symmetry in the covariance

    % Here are the partial lags
    ydim=[0:NyNx(1)-1]';
    xdim=[0:NyNx(2)-1] ;
    
    % Here is the Matern spatial covariance on the quarter distance grid, see above
    [Cyy,tyy]=spatmat(ydim,xdim,th,params,xver);

    % Exploit the symmetry just a tad, which allows us to work with smaller matrices
    q1=fft2(Cyy);
    q4=q1+[q1(:,1) fliplr(q1(:,2:end))];

    % Here is the blurred covariance on the 'complete' grid
    Hh=fftshift(2*real(q4-repmat(fft(Cyy(:,1)),1,NyNx(2)))...
	        -repmat(2*real(fft(Cyy(1,1:end))),NyNx(1),1)...
	        +Cyy(1,1));
    % If you ever wanted tyy/Cyy to come out you'll need to unquarter it
    tyy=[fliplr(tyy(:,2:end)) tyy]; tyy=[flipud(tyy(2:end,:)) ; tyy];
    tyy=[zeros(size(tyy,1)+1,1) [zeros(1,size(tyy,2)) ; tyy]];
    Cyy=[fliplr(Cyy(:,2:end)) Cyy]; Cyy=[flipud(Cyy(2:end,:)) ; Cyy];
    Cyy=[zeros(size(Cyy,1)+1,1) [zeros(1,size(Cyy,2)) ; Cyy]];
end

% Normalize and vectorize
Sbar=Hh(:)*prod(dydx)/(2*pi)^2;

% Check Hermiticity of the results
if xver==1
    blurcheck(Sbar,params)
    hermcheck(tyy)
    hermcheck(Cyy)
end

% Produce the unwrapped wavenumbers if you've requested them to be output
if nargout>1
  k=knums(params);
  k=k(:);
else
  k=[];
end

% Optional output
varns={Sbar,k,tyy,Cyy};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cyy,t]=spatmat(ydim,xdim,th,params,xver)
% [Cyy,t]=spatmat(ydim,xdim,th,params,xver)
%
% Returns the modified spatial covariance whose Fourier transform is the blurred
% spectrum after spatial data tapering, i.e. the expected periodogram. No use to
% return the autocovariance of the applied spatial taper which is an essential
% part of this operation...  never need it explicitly, used in SIMULOSL and
% LOGLIOSL via the intermediary of BLUROSY.

% Dimensions of the original grid
NyNx=params.NyNx;
dydx=params.dydx;

% Specify the spatial taper EXPLICITLY
if prod(size(params.taper))>1
    % Now compare with the other mechanisms
    % Completely general windows where 1 means you are taking a sample
    % See Arthur's note for more general windows, use IFF2/FFT2 you need, see
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/NonParametricEstimation/Periodogram.m
    % $MFILES/retired/QR?.m/map_*.m/whittle_*/expected_* etc
    Tx=double(params.taper);

    % If you are here with efs the taper is explicit AND not
    % symmetric, so must do something else
    if all([length(ydim) length(xdim)]==NyNx)
        % Produce the autocorrelation sequence eq. (12)
        t=zeros(size(Tx));
        % It's quite vital that these be colon ranges (faster) or (like
        % here) ROW index vectors... mixing rows/columns won't work
        for i=ydim(:)'+1
            for j=xdim(:)'+1
                % Vectorize? Check out XCORR2, that's also good
                t(i,j)=sum(sum(Tx(1:NyNx(1)-i+1,1:NyNx(2)-j+1).*(conj(Tx(i:end,j:end)))));
            end
        end
    else
        % This also obviates the extra test below, really this should be the
        % top way, but leave the explicit way for illustration
        t=xcorr2(Tx);
	% Add a row of zeros here
	t=[zeros(size(t,1)+1,1) [zeros(1,size(t,2)) ; t]];
    end
    % Now normalize the cross-correlations at the end
    t=t/sum(sum(Tx.^2));

    % Here too should use FFT where we can, see compute_kernels
    % internally and below, I would image that's just the same thing
    if xver==1 && all(size(t)==NyNx)
        t3=xcorr2(Tx); t3=t3/sum(sum(Tx.^2));
        if all(size(t)==NyNx)
            % Then the doubling inside this block needs to be undone
            t3=t3(NyNx(1):end,NyNx(2):end);
        end
        diferm(t,t3);
    end
else
    % Specify the spatial taper IMPLICITLY, taper is just a single number, could
    % be 0 or 1 both telling us "not" spatially tapered which is of course in
    % essence the same as a "unit" spatial taper taper operation. The triangle
    % functions are the normalized autocorrelations of the unit window functions,
    % i.e. c_{g,n}(u) of (12)-(13) in Guillaumin (2022), doi: 10.1111/rssb.12539
    triy=1-abs(ydim)/NyNx(1);
    trix=1-abs(xdim)/NyNx(2);
    % Here is the gridded triangle for this case
    t=bsxfun(@times,triy,trix);

    if xver==1
        % Do form the taper explicitly after all, normalize ahead of time
        Tx=ones(NyNx)/sqrt(prod(NyNx));
        % Need to cut one off Arthur says, possibly need to
        % re-re-visit these even/odd comparisons in BLUROS, if it ever
        % gets to that point; currently the comparison is favorable
        t2=fftshift(ifft2(abs(fft2(Tx,2*size(Tx,1)-1,2*size(Tx,2)-1)).^2));
        % Fix the rim by adding zeroes top and left
        t2=[zeros(size(t2,1)+1,1) [zeros(1,size(t2,2)) ; t2]];
        % Check the difference between these two implementations,
        % all checked for even/odd/method combinations on 2/24/2023
        if all(size(t)==NyNx)
            % Then the doubling inside this block needs to be undone
            t2=t2(NyNx(1)+1:end,NyNx(2)+1:end);
        end
        diferm(t,t2);
    end
end

% Here is the distance grid, whose size depends on the input
y=sqrt(bsxfun(@plus,[ydim*dydx(1)].^2,[xdim*dydx(2)].^2));

% The modified spatial covariance
Cyy=maternosy(y,th).*t;

% Remind me: Arthur: diag(U^T * Cyy * U) where U is the DFMTX is the
% expected periodogram, the diagonal variance, but to get the variance of
% the gradient we need the whole thing, the covariance, (U^T * Cyy * U) of
% course in two dimensions, so properly arranged.
