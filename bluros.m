function [Sbar,k,Fejk]=bluros(S,params,xver)
% [Sbar,k,Fejk]=BLUROS(S,params,xver)
%
% Blurring of a spectral matrix with the periodogram of a spatial
% windowing function (for now: the boxcar), in the approximate,
% discretized convolutional manner. The wavenumber-dependent input is
% given on a grid that was specified to be an integer refinement from
% an original that remains the target. The result is obtained by
% interpolation or subsampling to the original grid. This function is
% designed to be called only by MATERNOSP.
%
% INPUT:
%
% S       The spectral matrix with its wavenumbers unwrapped, on a
%         wavenumber grid that was REFINED from the original target
%         For univariate spectra, it's a column vector. 
%         For bivariate spectra, it's [S_XX(k(:)) SXY(k(:)) SYY(k(:))]
% params  Parameters of this experiment, the ones that are needed are:
%         dydx  sampling interval in the y and x directions [m m]
%         NyNx  number of samples in the y and x directions
%         blurs 0 Don't blur likelihood using the Fejer window
%               1 With no refinement, this is like not blurring at all!
%               N Blur likelihood using the Fejer window [default: N=2]
%              -1 Blur likelihood using the exact BLUROSY procedure
% xver    1 Extra verification, among other ways, via BLURCHECK
%         0 No checking at all
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, interpolated to the ORIGINAL requested
%         dimension as identified by 'params' in the input
% k       The wavenumber matrix (the norm of the wave vectors), unwrapped
% Fejk    The kernel, for now, only the Fejer kernel
%
% SEE ALSO:
%
% MATERNOSP, SIMULOSL, BLUROSY, BLUROSY_DEMO
%
% Last modified by fjsimons-at-alum.mit.edu, 11/26/2018

% If you tried running it as a stand-alone with the wrong parameters
if params.blurs<0
  error('You should be running BLUROSY, not BLUROS!')
end

% Set defaults
defval('xver',1)

% Target dimensions, i.e. the ORIGINAL ones
NyNx=params.NyNx;
dydx=params.dydx;

% The blurring REFINEMENT factor
blurs=params.blurs;

% The unblurred-property ORIGINAL wavenumber grid
[k,dci,~,kx,ky]=knums(params);
% The blurred-property REFINED wavenumber grid
[~,kzero,~,kx2,ky2]=knums(params,1);
% Find out what you've actually been given, do not jump to
% conclusions as KNUMS could introduce dimensional subtlety; at one
% point, we played with parity preservation, since abandoned
NyNx2=[length(ky2) length(kx2)];

if xver==1
  % Other checks are inside KNUM2 called from KNUMS called here
  diferm(k(dci(1),dci(2)))
  % Is the zero-wavenumber indeed, well, zero?
  diferm(k(kzero))
end

% This is the periodogram of the boxcar, for the REFINED grid as
% interpolated from the ORIGINAL grid (or we'd just hit the nodes).
% See http://blinkdagger.com/matlab/matlab-fft-and-zero-padding/
% In our subsequent formalism we force the fft/ifft pair to be unitary.
% If blurs were to be 1, we would get a single 1 in the center.
% Note that the kernel values are very different depending on even/odd dimensions
Fejk=fftshift(abs(fft2(repmat(1/sqrt(prod(NyNx))/sqrt(prod(NyNx2)),NyNx),NyNx2(1),NyNx2(2))).^2);

if xver==1
  % Make sure the kernel is unitary and norm-preserving
  difer(sum(Fejk(:))-1,[],[],NaN)
  % Check Hermiticity of the Fejer kernel, this NEVER fails
  hermcheck(Fejk)
end

% Check that there is no roundoff going on in the interpolation
% This case might apply when we start from an odd number and double
if [kx(end)-kx2(end)]>0; kx2(end)=kx(end); end
if [ky(end)-ky2(end)]>0; ky2(end)=ky(end); end

% This case might apply when we start from an even number and double
if [kx(1)-kx2(1)]<0; kx2(1)=kx(1); end
if [ky(1)-ky2(1)]<0; ky2(1)=ky(1); end

% Prefill to original dimensions
Sbar=nan(prod(NyNx),size(S,2));

% Remember that if S is a matrix, its three column vectors are
% elements of their joint spectral matrix, SXX(:), SXY(:), SYY(:)
for in=1:size(S,2)
  % Perform the convolution, isolate the center part where the
  % center is to be understood as the location of the
  % zero-wavenumber, see KNUMS and KNUM2
  HhC=conv2(Fejk,reshape(S(:,in),NyNx2),'same');

  % You just need to supply the correct zero-wavenumber value,
  % while it gets progressively better for larger patches and
  % higher refinement factors, it remains different in the
  % exact and the approximate methods. We need to do this outside,
  % in MATERNOSP, since that's the only place we have access to the
  % Matern parameters that tell us EXACTLY what is going on.

  % CORRECTION: MAYBE IT'S RIGHT? STILL ARGUES WE MIGHT EXCLUDE
  % THE ZERO WAVENUMBER, RATHER REPLACE THE ZERO WAVENUMBER BY THE
  % CORRECTLY SCALED VERSION WHICH IS VIA THE SPATIAL AVERAGE OF
  % THE WINDOW FUNCTION... 
  
  % disp(sprintf('\nINTERPOLATE:')) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tic
  disp(sprintf('Interpolating %ix%i to %ix%i',length(ky2),length(kx2),length(kx),length(ky)))
   % Interpolate, maybe later consider using griddedInterpolant
   HhI=interp2(kx2(:)',ky2(:),HhC,...
               kx(:)' ,ky(:) ,'linear');
    % The result, the convolution of two positive definite kernels,
    % cannot be negative! Happens under 'cubic' or 'spline'...
   if any(HhI<0); error('Convolved variance kernel contains negatives'); end
  d=toc;
  
  % That is it for now, unless the below changes it
  Hh=HhI;

  % disp(sprintf('\nSUBSAMPLE:')) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % There may be a simpler way in which the new grid is a superset of the old
  sx=1+mod(NyNx(2),2);
  sy=1+mod(NyNx(1),2);

  if round(blurs)==blurs && ...
          [sum(abs(kx2(sx:blurs:end)-kx))+sum(abs(ky2(sy:blurs:end)-ky))]==0
    tic
     HhS=HhC(sx:blurs:end,sy:blurs:end);
    e=toc;

    if e<d
      disp(sprintf('\nSUBSAMPLING FASTER THAN INTERPOLATION by %i%%',round((1-e/d)*100)))
    end

    % Check the difference
    disp('Relative difference between interpolated and subsampled versions')
    difer([HhI-HhS]/norm(HhS))

    % Isn't subsampling always to be preferred over interpolation?
    Hh=HhS;
  end

  % Need to consider that we can subsample in even more cases

  % We adjust for the bits that aren't already doubly present, 
  % watching where they are coming from, i.e. at the Nyquist
  if ~mod(NyNx(1),2)
    Hh(1,:)=Hh(1,:)*2;
  else
       Hh(end,:)=Hh(end,:)*2; 
  end
  if ~mod(NyNx(2),2)
    Hh(:,1)=Hh(:,1)*2; 
  else
     Hh(:,end)=Hh(:,end)*2; 
  end
  % Unwrap
  Sbar(:,in)=Hh(:);

  % Make sure the zero wavenumber gets the correct value, the zero-wavenumber
  if xver==1
    % We need the indicial zero-wavenumber location for the REFINED grid
    kzx2=floor(NyNx2(2)/2)+1; kzy2=floor(NyNx2(1)/2)+1;
    % We need the running zero-wavenumber location for the REFINED grid
    kz2=[kzx2-1]*NyNx2(1)+kzy2;
    % I am not sure yet which one of these two things it should be yet
    %    disp('Checking center portion')
    %diferm(Sbar(kzero,in),S(kz2,in)*Fejk(kzy2,kzx2)); 
    % Playing with this in order to find out as part of BLUROSY_DEMO
    disp(sprintf('\nArea factor %g\n',1/(1/sqrt(prod(NyNx))/sqrt(prod(NyNx2)))))
  end
end

if xver==1
  % Check that no extrapolation was demanded, effectively
  % but know that griddedInterpolant would have EXTRApolated fine
  difer(sum(isnan(Sbar(:))),[],2,NaN)
  % Check Hermiticity and positive-definiteness
  blurcheck(Sbar,params)
end

% Produce the unwrapped wavenumbers if you've request them to be output
if nargout>1
  k=k(:);
end

% disp(sprintf('BLUROS %i %i %i',blurs,NyNx2(1),NyNx2(2)));

% FURTHER CONSIDERATIONS PRECEDING THE MOVE TO BLUROSY
% Use the Claerbout helix? No gains, as once discussed with Sergey Fomel.
% Use CONVMTX2? Might need more memory
% User FFTW?
% Save and load the window?
% Allow for windows different from the boxcar? 
% Limit to halfplane? Now that's effectively done by BLUROSY

