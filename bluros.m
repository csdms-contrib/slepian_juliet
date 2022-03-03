function varargout=bluros(S,params,xver)
% [Sbar,k,Fejk]=BLUROS(S,params,xver)
%
% Blurring of a spectral matrix with the periodogram of a spatial windowing
% function (for now: the boxcar), in the approximate, discretized
% convolutional manner. The wavenumber-dependent input is given on a grid
% that was specified to be an integer refinement from an original that
% remains the target. The result is obtained by subsampling (or
% interpolation) to the original grid. This function is designed to be
% called only by MATERNOSP or BLUROSY_DEMO, where it is exhaustively tested.
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
% MATERNOSP, SIMULOSL, BLUROSY
%
% EXAMPLE:
%
% BLUROS('demo1',pp,bb) compares against BLUROSY where
%                pp     A square matrix size (one number)
%                   bb  A MATERNOSP blurring densification (one number)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2022

if ~isstr(S)
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
  [k,dci,dcn,kx,ky]=knums(params);
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
  % Note that the kernel values are very different depending on even/odd
  % dimensions. But it's not like it ever needs to be IFFTSHIFT.
  Fejk=fftshift(abs(fft2(repmat(1/sqrt(prod(NyNx))/sqrt(prod(NyNx2)),NyNx),NyNx2(1),NyNx2(2))).^2);

  if xver==1
    % Make sure the kernel is unitary and norm-preserving
    difer(sum(Fejk(:))-1,[],[],NaN)
    % Check Hermiticity of the Fejer kernel, this NEVER fails
    hermcheck(Fejk)
    % Check Hermiticity of the input spectrum, this NEVER fails
    for in=1:size(S,2)
      hermcheck(reshape(S(:,in),NyNx2))
    end
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

    % If F and S don't commute?
    % When the dimensions are even this fails to be Hermitian
    % but we're not concerned as the result is - we do see roundoff
    % disp(sprintf('Checking HhC for dimensions %ix%i',size(HhC)))
    % hermcheck(HhC)
    % disp(sprintf('%s\n','No message if passed'))
    % log10(HhC) % Seems to be not perfect but at least appropriate
    
    % This works out how the blurred grid is a superset of the blurred grid!
    sx=1+mod(NyNx(2),2)*floor(blurs/2);
    sy=1+mod(NyNx(1),2)*floor(blurs/2);
    
    if round(blurs)==blurs &&  ...
          [sum(abs(kx2(sx:blurs:end)-kx))+sum(abs(ky2(sy:blurs:end)-ky))]<eps
      % disp(sprintf('\nSUBSAMPLE:')) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp(sprintf('Subsampling %ix%i to %ix%i',length(ky2),length(kx2),length(kx),length(ky)))
      HhS=HhC(sx:blurs:end,sy:blurs:end);
      % Isn't subsampling always to be preferred over interpolation?
      Hh=HhS;
    else
      % disp(sprintf('\nINTERPOLATE:')) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp(sprintf('Interpolating %ix%i to %ix%i',length(ky2),length(kx2),length(kx),length(ky)))
      % Interpolate, maybe later consider using griddedInterpolant
      HhI=interp2(kx2(:)',ky2(:),HhC,...
                  kx(:)' ,ky(:) ,'nearest');
      % The result, the convolution of two positive definite kernels,
      % cannot be negative! Happens under 'cubic' or 'spline'...
      if any(HhI<0); error('Convolved variance kernel contains negatives'); end
      
      % That is it for now, unless the below changes it
      Hh=HhI;

      % Keep the comparison for while you're developing
      % disp('Relative difference between interpolated and subsampled versions')
      % diferm([HhI-HhS]/norm(HhS))
    end

    % We adjust for the bits that aren't already doubly present, 
    % watching where they are coming from, i.e. at the Nyquist
    % which is only at the top and left for even dimensions
    % And watch the zero? See RANDGPN? Not quite there yet, are we
    % but in the comparison that's definitely the ticket
    % As in - double all the points that don't have symmetric
    % counterparts, that's all, isn't it. But not the zero lines since
    % these already existed in the original. 
    d=0;
    if ~mod(NyNx(1),2)
      Hh(1,:)=Hh(1,:)*2;
      d=d+1;
    end
    if ~mod(NyNx(2),2)
      Hh(:,1)=Hh(:,1)*2;
      d=d+1;
    end
    % Nooooo ; if d==2; Hh(1,1)=Hh(1,1)/2; end

    % No further adjustments necessary (though I tried corner points etc)
    % and frankly, got a little tired of it. Check details of (16,4) and
    % (16,5) to see perhaps what's going on
    
    % Should we make Hh absolutely super-duper Hermitian as a
    % precautionary measure?

    % So if the input is e/o and the refinement o, no parity change, all
    % good, but if the input is e/o and the refinement e, parity changes for
    % the odd and then going back from even to odd messes things up slightly.
    
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
      % diferm(Sbar(kzero,in),S(kz2,in)*Fejk(kzy2,kzx2));
      % Sbar(kzero,in)=S(kz2,in)*Fejk(kzy2,kzx2);
      % Playing with this in order to find out as part of BLUROSY_DEMO
      % disp(sprintf('\nArea factor %g\n',1/(1/sqrt(prod(NyNx))/sqrt(prod(NyNx2)))))
    end
  end

  if xver==1
    % Check that no extrapolation was demanded, effectively
    % but know that griddedInterpolant would have EXTRApolated fine
    difer(sum(isnan(Sbar(:))),[],2,NaN)
    % Check Hermiticity and positive-definiteness - later loop for multi-D
    disp(sprintf('Checking Sbar for dimensions %ix%i',params.NyNx))
    blurcheck(Sbar,params)
    disp(sprintf('%s\n','No message if passed'))
    
    % Remember if it's even we disregard the Nyquist for symmetry
    A=Hh(1+~mod(size(Hh,1),2):end,1+~mod(size(Hh,2),2):end);

    % Some on-screen eyeballing
    % conj(fliplr(flipud(A)))-A
    % conj(fliplr(flipud(A)))./A
    
    % At some point this was a testing tool
    % figur(3)
    % subplot(211)
    % %  imagesc(log10(abs(conj(fliplr(flipud(A)))-A))); axis image; colorbar
    % imagesc(conj(fliplr(flipud(A)))-A); axis image; colorbar
    % title('Additive non-symmetry of the blurred spectral matrix')

    % subplot(212)
    % imagesc(conj(fliplr(flipud(A)))./A); axis image; colorbar
    % title('Multiplicative non-symmetry of the blurred spectral matrix')
  end

  % Variable output
  varns={Sbar,k(:),Fejk};
  varargout=varns(1:nargout);

  % disp(sprintf('BLUROS %i %i %i',blurs,NyNx2(1),NyNx2(2)));

elseif strcmp(S,'demo1')
  % Just call the formerly separate function, not subfunction, with the
  % appropriate parameters passed on through
  % Remember the parameter names are no longer what they seem, the first
  % one is pp and the second one bb
  blurosy_demo(params,xver)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fejk=blurosy_demo(pp,bb)
% Fejk=BLUROSY_demo(pp,bb)
%
% This function compares the various ways of blurring, make sure you
% run it in a variety of combinations in both even and odd numbers.
%
% INPUT:
%
% pp       A square matrix size (one number, e.g. 32, or randi(100))
% bb       A MATERNOSP blurring densitfication (e.g. 3, or randi(10))
%
% OUTPUT:
%
% Fejk    The kernel, for now, only the Fejer kernel
%
% Last modified by fjsimons@alum.mit.edu, 03/3/2022

% Some random input values, but the product must be high!
defval('pp',randi(100))
defval('bb',randi(7))

% Start with an even p.NyNx, make the odd later
pp=2*round(pp/2);

% Supply the grid size and spacing and the Matern covariance parameters
p.NyNx=[pp pp];
p.dydx=1e3*[1 1]; 
th=1e6*[1 0.0000025 0.02]; 

% Perform the blurring of the spectral matrix using different techniques
% Convolutionally blurred with the specified refinement parameter
p.blurs= bb; S2=maternosp(th,p,1); 
 % Exactly blurred, slow (do not use MATERNOSP; it defaults to 'efs')
p.blurs=-1; S3=blurosy(th,p,1,'ef');
% Exactly blurred, fast (do not use MATERNOSP; it defaults to 'efs')
p.blurs=-1; S4=blurosy(th,p,1,'efs');

% Make the plots for visual inspection! 
figure(1)
clf
[ah,ha,H]=krijetem(subnum(2,3)); 

% Inspect the blurred spectra on a decibel scale (axis adjusted later)
axes(ah(1))
imagesc(reshape(decibel(S2),p.NyNx)); 
t(1)=title(sprintf('[%i %i] BLUROS-%i   [S2]',p.NyNx,bb));
axes(ah(2))
imagesc(reshape(decibel(S3),p.NyNx)); 
t(2)=title(sprintf('[%i %i] BLUROSY-ef  [S3]',p.NyNx));
axes(ah(3))
imagesc(reshape(decibel(S4),p.NyNx)); 
t(3)=title(sprintf('[%i %i] BLUROSY-efs [S4]',p.NyNx));

% Move that title
movev(t,-pp/20)

% Plot the zero-wavenumber as a cross, to keep track
[kor,dci,dcn,kx,ky]=knums(p); for in=1:3; axes(ah(in));
hold on; plot(dci(1),dci(2),'k+'); axis image ; hold off;
set(ah(in),'xtick',unique([1 dci(2) p.NyNx(2)]),...
	   'ytick',unique([1 dci(1) p.NyNx(1)])); end

% Max of the blurred spectrum needs to be at the zero wavenumber,
% or else it will complain
kz=[dci(2)-1]*p.NyNx(1)+dci(1);
[m2,i2]=max(S2); diferm(i2-kz)
[m3,i3]=max(S3); diferm(i3-kz)
[m4,i4]=max(S4); diferm(i4-kz)

%disp(sprintf('\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Moving into the plus-one territory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'))

% Whatever you just did, add one to the dimension to hit the
% odd/even pairs
p.NyNx=p.NyNx+1;
% Convolutionally blurred
p.blurs= bb; S21=maternosp(th,p,1); 
% Exact blurred slow
p.blurs=-1; S31=blurosy(th,p,1,'ef'); 
% Exact blurred fast
p.blurs=-1; S41=blurosy(th,p,1,'efs'); 

axes(ah(4))
imagesc(reshape(decibel(S21),p.NyNx)); 
t2(1)=title(sprintf('[%i %i] BLUROS-%i   [S21]',p.NyNx,bb));
axes(ah(5))
imagesc(reshape(decibel(S31),p.NyNx))
t2(2)=title(sprintf('[%i %i] BLUROSY-ef  [S31]',p.NyNx));
axes(ah(6))
imagesc(reshape(decibel(S41),p.NyNx)); 
t2(3)=title(sprintf('[%i %i] BLUROSY-efs [S41]',p.NyNx));

% Move that title
movev(t2,-pp/20)

% Wavenumber crosses
[kor,dci,dcn,kx,ky]=knums(p); for in=4:6; axes(ah(in));
hold on; plot(dci(1),dci(2),'k+'); axis image ; hold off; 
set(ah(in),'xtick',unique([1 dci(2) p.NyNx(2)]),...
           'ytick',unique([1 dci(1) p.NyNx(1)])); end

% Max needs to be at the zero wavenumber!
kz2=[dci(1)-1]*p.NyNx(1)+dci(2);
[m21,i21]=max(S21); diferm(i21-kz2)
[m31,i31]=max(S31); diferm(i31-kz2)
[m41,i41]=max(S41); diferm(i41-kz2)

% Just tell us what the maxima are, so we can fix down the line
disp('Rounded maxima of every panel')
disp(sprintf('%i: %i %i %i\n%i: %i %i %i',...
	     round([pp m2 m3 m4 pp+1 m21 m31 m41])))

% Total amounts of spectral energy, somehow?
% disp(sprintf('%i: %i %i %i\n%i: %i %i %i',...
%	     round([pp sum(S2(:)) sum(S3(:)) sum(S4(:)) pp+1 sum(S21(:)) ...
%                    sum(S31(:)) sum(S41(:))])))

% Align plots
longticks(ah); serre(H',1/2,'down')
for in=1:6; axes(ah(in)); caxis([-75 0]); end

% Print it baby
figdisp([],sprintf('1_%3.3i_%2.2i',pp,bb),[],2)

% Second plot that compares the ratios of the computed quantities
figure(2)
clf
[ah2,ha2,H2]=krijetem(subnum(2,2));
axes(ah2(1))
plot(S2./S3,'.'); 
t3(1)=title(sprintf('[%i %i] BLUROS-%i S2/S3',[pp pp],bb));
hold on; plot(i2,S2(i2)/S3(i3),'o'); hold off
set(ah2(1),'xtick',unique([1 kz pp^2]),'xgrid','on')
axes(ah2(2))
plot(S3./S4,'.'); 
t3(2)=title(sprintf('[%i %i] S3/S4',[pp pp]));
hold on; plot(i3,S3(i3)/S4(i4),'o'); hold off
set(ah2(2),'xtick',unique([1 kz pp^2]),'xgrid','on')
axes(ah2(3))
plot(S21./S31,'.'); 
t3(3)=title(sprintf('[%i %i] BLUROS-%i S21/S31',[pp pp]+1,bb));
hold on; plot(i21,S21(i21)/S31(i31),'o'); hold off
set(ah2(3),'xtick',[1 kz2 [pp+1]^2],'xgrid','on')
axes(ah2(4))
plot(S31./S41,'.'); 
t3(4)=title(sprintf('[%i %i] S31/S41',[pp pp]+1));
hold on; plot(i31,S31(i31)/S41(i41),'o'); hold off
set(ah2(4),'xtick',[1 kz2 [pp+1]^2],'xgrid','on')

% Align plots
longticks(ah2); 

for index=1:length(ah2)
  axes(ah2(index))
  axis tight
  yl=ylim; xl=xlim;
  hold on; plot(xl,[1 1],':')
  % Better diagnosing
  set(ah2(index),'ylim',[min(0.75,min(yl)*0.9)-[min(yl)<range(yl)/20]*range(yl)/20 max(yl)*1.1])
  set(ah2(index),'ylim',[0 2])
end

% Print it baby
figdisp([],sprintf('2_%3.3i_%2.2i',pp,bb),[],2)

% A CONSIDERATION AT SOME POINT BUT NOW LARGELY FORGOTTEN

% Do we need to supply the correct zero-wavenumber value? While it gets
% progressively better for larger patches and higher refinement factors, it
% remains different in the exact and the approximate methods. We need to do
% this outside, in MATERNOSP, since that's the only place we have access to
% the Matern parameters that tell us EXACTLY what is going on. Maybe replace
% the zero wavenumber by the correctly scaled version which is via the
% spatial average of the window function... 

% FURTHER CONSIDERATIONS PRECEDING THE ULTIMATE MOVE TO BLUROSY
%
% Use the Claerbout helix? No gains, as once discussed with Sergey Fomel.
% Use CONVMTX2? Might need more memory...
% User FFTW? ...
% Save and load the window? ...
% Allow for windows different from the boxcar? ...
% Limit to halfplane? Now that's effectively done by BLUROSY ...
