function varargout=blurosy(th,params,xver,method,tsto)
% [Sbar,k,tyy,Cyy]=blurosy(th,params,xver,method,tsto)
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
% tsto    An extra parameter slot to pass onto demo2    
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
% BLUROSY('demo1',pp,bb) compares against BLUROS where
%                 pp     A square matrix size (one number)
%                    bb  A MATERNOSP blurring densification (one number)
%
% BLUROSY('demo2',pp,nn,mm,tt) Boxcar/Ukraine/France example
%                 pp           A matrix size (two numbers) [defaulted]
%                    nn        A number of iterations to average over [defaulted] 
%                       mm     A method, either 'ef' or 'efs' [defaulted]
%                          tt  1 Boxcar 2 France 3 Ukraine
%
% BLUROSY('demo3') % should produce no output
%
% Last modified by arthur.guillaumin.14-at-ucl.ac.uk, 10/15/2017
% Last modified by fjsimons-at-alum.mit.edu, 09/22/2023

if ~isstr(th)
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
elseif strcmp(th,'demo1')
    % This case is handled by BLUROS
    % Remember the prompt parameter names are no longer what they seem, the first
    % one is matrix size and the second one blurring refinement density den
    bluros('demo1',params,xver)
elseif strcmp(th,'demo2')
    % Now here we will test that the blurred spectrogram is the average
    % periodogram of the periodogram of the sample generated by the spectral
    % density. We welcome variability around it, as this is in line with
    % the uncorrelated-wavenumber approximation of the debiased Whittle
    % approach, but we do fear systematic offsets.

    % Size of the patch
    defval('params',[33 33])
    % Number of iterations    
    defval('xver',100)
    % Method of computation
    defval('method','efs')
    % Type of taper, e.g. boxcar, France, Ukraine
    defval('tsto',1)
    
    % Some combinations - SIMULOSL
    p.NyNx=[params(1) params(2)];
    p.dydx=1e3*[1 1];
    % Some Matern paramters
    th=1e6*[1 0.0000025 0.001]; 

    % Compare the average periodogram with the blurred spectral density

    % What kind of a test are we running? Boxcar, or France/Ukraine?
    switch tsto
      case 1
        % Boxcar is 0 or 1 or ones
        p.taper=1;
        % The mask parameter is irrelevant, though it might be the anti-taper
      case 2
        % Here's another one
        p.mask='france';
      case 3
        % Here's another one
        p.mask='ukraine';
    end
    switch tsto
      case {2,3}
        % Generate "mask" only, never mind what the data will be
        [~,~,I]=maskit(rand(p.NyNx),p);
        % Keep the mask as a taper or the taper as mask, for illustration only
        p.taper=~I;
    end
    
    % Calculate expected periodogram, i.e. the appropriately blurred likelihood
    % Use the exact method via BLUROSY, force p.blurs=-1
    p.blurs=-1;
    % Just do the xver=1 explicitly here
    Sbar=blurosy(th,p,1,method);
    % Then for what comes next, to simulate data using SGP, force p.blurs=Inf
    p.blurs=Inf;

    % First figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    clf
    [ah,ha,H]=krijetem(subnum(2,3));
    
    % The spatial domain taper
    axes(ah(4))
    if length(p.taper)==1; p.taper=ones(p.NyNx); end
    imagesc(p.taper); axis image
    t(4)=title('taper');

    % The expectation of the periodogram
    axes(ah(3))
    imagesc(log10(v2s(Sbar,p))); axis image
    t(3)=title(sprintf('%s [%g %g %g] | expectation',...
                       '\theta =',th./[1 1 sqrt(prod(p.dydx))]))
    
    % We be collecting the average periodogram to compare with the expected periodogram
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
    % See EGGERS6/8 for cleanup
    subplot(211)
    hist(2*Sb./Sbar,1.5*log2(round(prod(p.NyNx))+1))
    shrink(gca,2,1)
    
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
        H.LineWidth=1;
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
elseif strcmp(th,'demo3')
    % Simulate some random data with default p.blurs=Inf and p.taper=0
    [H,th,p]=simulosl;
    % Reset to implicit treatment of unit taper
    p.blurs=-1; p.taper=1;
    % This function produces the blurred spectral densities, different methods
    [Sbar1,k1,tyy1,Cyy1]=blurosy(th,p,1,'ef');
    [Sbar2,k2,tyy2,Cyy2]=blurosy(th,p,1,'efs');
    % Reset to explicit treatment of unit taper
    p.taper=ones(p.NyNx);
    [Sbar3,k3,tyy3,Cyy3]=blurosy(th,p,1,'ef');
    [Sbar4,k4,tyy4,Cyy4]=blurosy(th,p,1,'efs');
    % All of these should be virtually identical
    % Note the Sbar are large so the tolerance is high 
    diferm(Sbar1,Sbar2,-2); diferm(tyy1,tyy2); diferm(Cyy1,Cyy2);
    diferm(Sbar1,Sbar3,-2); diferm(tyy1,tyy3); diferm(Cyy1,Cyy3);
    diferm(Sbar1,Sbar4,-2); diferm(tyy1,tyy4); diferm(Cyy1,Cyy4);
    diferm(Sbar2,Sbar3,-2); diferm(tyy2,tyy4); diferm(Cyy2,Cyy3);
    diferm(Sbar2,Sbar4,-2); diferm(tyy2,tyy4); diferm(Cyy2,Cyy4);
    diferm(Sbar3,Sbar4,-2); diferm(tyy3,tyy4); diferm(Cyy3,Cyy4);
end

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
