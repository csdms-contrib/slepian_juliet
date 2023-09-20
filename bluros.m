function varargout=bluros(S,params,xver,tsto)
% [Sbar,k,c2,Fejk]=BLUROS(S,params,xver,tsto)
%
% Blurring of a spectral matrix with the periodogram of a spatial
% windowing function (for now: the boxcar, hence the Fejer kernel), in
% the approximate, discretized convolutional manner. The
% wavenumber-dependent input is given on a grid that was specified to
% be an integer refinement from an original that remains the
% target. The result is obtained by subsampling (or interpolation) to
% the original grid. This function is called only by MATERNOSP or by
% itself for demoing and testing, which has been exhaustive. 
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
%         blurs 0 No wavenumber blurring
%               1 No wavenumber blurring, effectively
%               N Fejer convolutional blurring on an N-times refined grid
%              -1 Fejer multiplicative BLUROSY using exact procedure, error
%             Inf Simulate using SGP invariant embedding, error (N/A!)
%         taper  0 there is no taper near of far
%                1 it's a unit taper, implicitly
%                OR an appropriately sized taper with proper values 1 or 0
%                   not supported here yet, for that use BLUROSY   
% xver    1 Extra verification, among other ways, via BLURCHECK
%         0 No checking at all
% tsto    An extra parameter slot to pass onto demo2    
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, interpolated to the ORIGINAL requested
%         dimension as identified by 'params' in the input
% k       The wavenumber matrix (the norm of the wave vectors), unwrapped
% c2      The coupling kernel that I think is (A58) in Olhede & Simons (2013)
% Fejk    The Fejer kernel
%
% SEE ALSO:
%
% MATERNOSP, SIMULOSL, BLUROSY
%
% EXAMPLE:
%
% BLUROS('demo1',pp,bb) compares against BLUROSY where
%                pp     A square-matrix size (one number)
%                   bb  A MATERNOSP blurring densification (one number)
%
% BLUROS('demo2',pp,nn,tt) Boxcar/Ukraine/France example
%                pp        A matrix size (two numbers)
%                   nn     A number of iterations over which will be averaged
%                      tt  1 Boxcar 2 France 3 Ukraine
%
% BLUROS('demo3',33,3)
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/17/2023

if ~isstr(S)
  % disp('Ich bin so lang nicht bei dir g''west')
  % If you tried running it as a stand-alone with the wrong parameters
  if params.blurs<0 | isinf(params.blurs)
    error('You should be calling BLUROSY, not BLUROS!')
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

  defval('c2',[])
  if nargout>2
    % 1-D alternatives and using L'Hopital
    kyy=ky2/2*params.dydx(1); kxx=kx2/2*params.dydx(2);
    Dy=sin(kyy*NyNx(1))./sin(kyy)/sqrt(NyNx(1))/sqrt(NyNx2(1));
    Dy(isnan(Dy))=sqrt(NyNx(1))/sqrt(NyNx2(1));
    Dx=sin(kxx*NyNx(2))./sin(kxx)/sqrt(NyNx(2))/sqrt(NyNx2(2));
    Dx(isnan(Dx))=sqrt(NyNx(2))/sqrt(NyNx2(2));
    % 2-D alternative for comparison with Fejk
    % a=repmat(Dy(:).^2,1,NyNx2(2));
    % b=repmat(Dx(:)'.^2,NyNx2(1),1);
    % This is much faster obviously
    Dirkk=Dy(:)*Dx(:)';
    % Check the difference
    diferm(Fejk-Dirkk.^2)
    % Now play with the convolution of the product of two Dirichlet kernels?
    DD=conv2(Dirkk,Dirkk,'same');
    keyboard
    % D2=DD.^2;
    % Maybe THAT is what we should use on the inside of the Fisher
    % approximation, first subsample dd
    sx=1+mod(NyNx(2),2)*floor(blurs/2);
    sy=1+mod(NyNx(1),2)*floor(blurs/2);
    dd=DD(sx:blurs:end,sy:blurs:end);
    % All the other ones - clearly underexploiting symmetry
    sp=1:prod(NyNx);
    % Unwrapped k versus unwrapped k'
    for s=1:prod(NyNx)
      % Calculate the offsets from the center
      [~,di,dj]=sspdist(s,sp,NyNx);
      % Then fill them all in
      % Decide what to do when the indices exceed what you have
      % available. Maybe periodic or mod somehow.
      noy=[dci(1)+di]>0 & [dci(1)+di]<=NyNx(1);
      nox=[dci(2)+dj]>0 & [dci(2)+dj]<=NyNx(2);
      they=dci(1)+di;
      thex=dci(2)+dj;
      con=noy & nox;
      % This is the coupling matrix...
      c2(s,sp(con))=dd(sub2ind(NyNx,they(con),thex(con))).^2;
      % Only constant diagonal when NyNx is odd?
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

  % disp(sprintf('BLUROS %i %i %i',blurs,NyNx2(1),NyNx2(2)));

elseif strcmp(S,'demo1')
  % Just call the formerly separate function, not subfunction, with the
  % appropriate parameters passed on through
  % Remember the parameter names are no longer what they seem, the first
  % one is pp and the second one bb
    blurosy_demo(params,xver)
    % You're not producing any output
    [Sbar,k,c2,Fejk]=deal([]);
elseif strcmp(S,'demo2')
  % Now here we will test that the blurred spectrogram is the average
  % periodogram of the periodogram of the sample generated by the spectral
  % density. We welcome variability around it, as this is in line with
  % the uncorrelated-wavenumber approximation of the debiased Whittle
  % approach, but we do fear systematic offsets.

  % Size of the patch
  defval('params',[33 33])
  % Number of iterations    
  defval('xver',100)
  % Type of taper
  defval('tsto',1)
    
  % Some combinations - SIMULOSL
  p.NyNx=[params(1) params(2)];
  p.dydx=1e3*[1 1]; 
  th=1e6*[1 0.0000025 0.001]; 

  % What kind of a test are we running? Boxcar, or France/Ukraine?
  switch tsto
    case 1
      % Compare the average periodogram with the blurred spectral density
      p.taper=1;
      p.mask=1;
    case 2
      % Here's another one
      p.mask='france';
      % Generate mask only
      [~,~,I]=maskit(rand(p.NyNx),p);
      p.taper=~I;
    case 3
      % Here's another one
      p.mask='ukraine';
      % Generate mask only
      [~,~,I]=maskit(rand(p.NyNx),p);
      p.taper=~I;
  end

  % So you're going to simulate using SGP and compare with BLUROSY
  p.blurs=-1;
  Sbar=blurosy(th,p);
  p.blurs=Inf;

  % First figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1)
  clf
  [ah,ha,H]=krijetem(subnum(2,3));
        
  % The spatial domain taper
  axes(ah(4))
  if length(p.taper)==1; p.taper=ones(p.NyNx); end
  imagesc(p.taper); axis image
  t(4)=title('mask');

  % The expectation of the periodogram
  axes(ah(3))
  imagesc(log10(v2s(Sbar,p))); axis image
  t(3)=title(sprintf('%s [%g %g %g] | expectation',...
                     '\theta =',th./[1 1 sqrt(prod(p.dydx))]))
  
  % muck it up to purposely break it
  %th(3)=th(3)*1.5
    
  Sbb=0;
  % The third input in the demo was the number of iterations, fake-called 'xver'
  for index=1:xver
      % Simulate sequentially and collect the average periodogram
      [Hx,th0,p,k,Hk,Sb,Lb,gane,miy]=simulosl(th,p);
      % One random one from the sequence will be shown
      randix=randi(index);
      if index==randix;
          axes(ah(1))
          imagefnan([1 1],p.NyNx([2 1]),v2s(Hx,p),[],halverange(Hx,75)); axis image ij
          t(1)=title(sprintf('Realization # %i',index));
          
          axes(ah(2))
          imagesc(log10(v2s(Sb,p))); axis image
          t(2)=title(sprintf('Realization # %i',index));
      end
      % Collect the average, watch the growth later
      Sbb=Sbb+Sb;
      % Mean of the evolving ratio over the wavenumbers
      m(index)=mean(Sbb/index./Sbar);
      % Mean of the evolving standard deviation of the mean ratio over the wavenumbers
      s(index)=std(Sbb/index./Sbar);
  end 
  % Mean ratio over the realizations
  Sbb=Sbb/xver;

  % Then compare with the thing coming out of BLUROSY
  axes(ah(6))
  imagesc(log10(v2s(Sbb,p))); axis image
  t(6)=title(sprintf('aver over %i realizations',xver));
      
  axes(ah(5))
  imagesc((v2s(Sbb./Sbar,p))); axis image
  t(5)=title(sprintf('(aver / expec), m %4.2f, s %4.2f',...
                     m(end),s(end)));
  try
      set(ah(5),'clim',m(end)+[-1 1]*2*s(end))
  end
  
  % Clean that sh*t up
  [k,dci,dcn,kx,ky]=knums(p);
  set(ah([1 4]),'xtick',unique([1 dci(2) p.NyNx(2)]),...
     'ytick',unique([1 dci(1) p.NyNx(1)]));
  set(ah([2 3 5 6]),'xtick',unique([1 dci(2) p.NyNx(2)]),...
     'ytick',unique([1 dci(1) p.NyNx(1)]),...
     'xticklabel',[-1 0 1],...
     'yticklabel',[-1 0 1]);
  longticks(ah)
  set(ah([3 6]),'YAxisLocation','right')
  % serre(H,1,'across')
  serre(H',0.5,'down')
  movev(t,-p.NyNx(1)/20)
  figdisp([],sprintf('demo_2_%3.3i_%3.3i_%3.3i',p.NyNx,xver),[],2)

  % Second figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(2)
  clf
  df=2;
  % See EGGERS8
  subplot(211)
  hist(2*Sb./Sbar,1.5*log2(round(prod(p.NyNx))+1))
  
  subplot(212)
  % This should be a straight line folks since the ratio is 1/2 chi^2_2
  h=qqplot(2*Sb./Sbar,makedist('gamma','a',df/2,'b',2)); axis equal;
  hx=get(h,'Xdata'); hx=hx{1};
  hy=get(h,'ydata'); hy=hy{1};
  qq0=plot([0 4*df],[0 4*df],'k'); hold on
  axis image
  plot(hx,hy,'LineS','none','Marker','o','MarkerF','r',...
       'MarkerE','r','MarkerS',2);
  longticks(gca); grid on
  xlim([0 4*df]); 
  ylim([0 4*df]); 
  %plot(2*Sbb./Sbar,1-exp(-2*Sbb./Sbar),'+'); axis tight normal
  %plot(2*Sbb./Sbar,chi2pdf(2*Sbb./Sbar,df),'+'); axis tight normal
  figdisp([],sprintf('demo_2_%3.3i_%3.3i_%3.3i_chi',p.NyNx,xver),[],2)

  % Second figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(3)
  clf
  H=errorbar(1:xver,m,-2*s,2*s);
  try
      set(getkids(H,1),'LineWidth',2)
      set(getkids(H,2),'Color',grey); 
  catch
      H.LineWidth=2;
      H.Color=grey;
  end
  longticks(gca,2)
  axis tight; grid on; xlabel('Sample size'); 
  ylabel('Ratio of average to expected periodogram')
  % set(gca,'xtick',[1 10:10:100]);
  shrink(gca,1,1.1)
  ylim([min(m-s*2)*1.1 max(m+2*s)*1.1])
  title(sprintf('%i x %i | %s = [%g %g %g]',p.NyNx,'\theta',...
                th./[1 1 sqrt(prod(p.dydx))]))
  figdisp([],sprintf('demo_2_%3.3i_%3.3i_%3.3i-%3.3i',p.NyNx,1,xver),[],2)

  % You're not producing any output
  [Sbar,k,c2,Fejk]=deal([]);
elseif strcmp(S,'demo3')
    % Some combinations - SIMULOSL
    % Input parsing for size
    p.NyNx=[params params];
    % Input parsing for blurring
    p.blurs=xver;
    
    p.dydx=1e3*[1 1]; 
    th=1e6*[1 0.0000025 0.001];

    % So this is what MATERNOSP actually does also, we just want the output of BLUROS
    % We're going to calculate the refinement convolution of two Dirichlet kernels
    [Sbar,k,c2]=bluros(maternos(knums(p,1),th),p);
    % Now inspect c2
    imagesc(log10(c2)); axis image
    colorbar
    % No output
    Fejk=NaN;
end

% Variable output
varns={Sbar,k(:),c2,Fejk};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function blurosy_demo(pp,bb)
% BLUROSY_demo(pp,bb)
%
% This function compares the various ways of blurring, make sure you
% run it in a variety of combinations in both even and odd numbers.
%
% INPUT:
%
% pp       A square matrix size (one number, e.g. 32, or randi(100))
% bb       A MATERNOSP blurring densitfication (e.g. 3, or randi(10))
%
% Last modified by fjsimons@alum.mit.edu, 01/24/2023

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

end

