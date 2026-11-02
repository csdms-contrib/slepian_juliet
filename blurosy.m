function varargout=blurosy(th,params,xver,method,tsto,dth)
% [SbarordSbardth,k,tyy,CyyordCyydth]=blurosy(th,params,xver,method,tsto,dth)
%
% Wavenumber blurring of a univariate Matern spectral density (or its spectral
% derivative) with the periodogram of a spatial taper. The result is the expected
% periodogram (or a term that enters the variance calculation).
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
% xver    1 extra verification via BLURCHECK and alternative computations
%         0 no checking at all
% method  'ef'  exact, efficient and fast [default]
%         'efs' exact, efficient and exploiting symmetry
%         'efd' exact, for diagonals calculation in JMATRIX method 3
%               which requires tsto but should be close to ef/efs is tsto=[0 0]
%               (except the result is fftshifted... )
%               (remaining question: does that solve BLUROS_DEMO lingering issue
%                   -> No, run BLUROSY DEMO1 to see comparison)
% tsto    An extra parameter slot to pass onto demo2
%         OR the 2 element list of offsets for the frequency indices, [m1, m2]
%            default assumption is [0, 0]; used for COVGAMMIOSL method 3
% dth     1, 2, or 3 specifies which element of th gets differentiated, to serve as
%         an input to the variance calculation, which requires a "blurred" version of
%         MAOSL, if you will, but there is no sense doing it there.
%
% OUTPUT:
%
% SbarordSbardth    The blurred spectral matrix or its derivative, with the unwrapped
%                   requested dimension as identified by the input 'params' (wrap with V2S)
% k                 The wavenumber matrix (the norm of the wave vectors), unwrapped
% tyy               The autocorrelation of the spatial taper... which you 
%                   may never need explicitly, used in SIMULOSL and LOGLIOSL
% CyyordCyydth      The equivalent modified Matern correlation, may never need it explicitly
%
% SEE ALSO:
%
% SIMULOSL, BLUROS, MATERNOSP, BLURCHECK, GETTAPER
%
% EXAMPLE:
%
% BLUROSY('demo1',pp,bb) compares against BLUROS where
%                 pp     A square matrix size (one number)
%                    bb  A MATERNOSP blurring densification (one number)
%
% BLUROSY('demo2',pp,nn,mm,tt,cc) Unit/Ukraine/France/Speckle/
                                  %Cosine/Rectangle/Structured example
%                 pp              A matrix size (two numbers) [defaulted]
%                    nn           A number of iterations to average over [defaulted]
%                                 OR: a date string (e.g. '14-Oct-2025') will produce
%                                   a second row of figures with the distribution of
%                                   the parameter (s_x^2) from a MLEOSL precalculated
%                                   ensemble that are either loaded or calculated live
%                                   (see also EGGERS6 for pedestrian tapers)
%                       mm        A method, either 'ef' or 'efs' [defaulted]
%                          tt     1 Unit    2 France    3 Ukraine    4 Speckle 
%                                 5 %Cosine 6 Rectangle 7 Structured 
%                                 8 SmoothFrance 9 <NAME>
%                             cc  Flag for providing a tapered dataset but
%                                 omitting the tapering in the analysis
%
% BLUROSY('demo3') % should produce no output
%
% Last modified by arthur.guillaumin.14-at-ucl.ac.uk, 10/15/2017
% Last modified by fjsimons-at-alum.mit.edu, 06/24/2025
% Last modified by olwalbert-at-princeton.edu, 09/08/2026

if ~isstr(th)
    % Defaults (avoiding DEFVAL to avoid costly EVALIN statements)
    if ~exist('xver','var') | isempty(xver); xver=1; end
    if ~exist('method','var') | isempty(method); method='ef'; end
    if ~exist('tsto','var'); tsto=[]; end
    if ~exist('dth','var'); dth=[]; end

    if params.blurs>=0 && ~isinf(params.blurs)
        error('Are you sure you should be running BLUROSY, not BLUROS?')
    end
    % Unlikely this will every trigger as exact blurring will be -1 
    if isinf(params.blurs) & isempty(dth)
        error('Are you sure you should be running BLUROSY, not MATERNOSY?')
    end

    % Target dimensions, the original ones
    NyNx=params.NyNx;
    dydx=params.dydx;

    switch method 
      case 'ef'
        % Generates a 2*NyNx double grid from which we subsample
        % Fully exact and not particularly fast, still much faster than BLUROS

        % Here are the full lags
        ydim=[-NyNx(1):NyNx(1)-1]';
        xdim=[-NyNx(2):NyNx(2)-1] ;

        % Here is the Matern spatial covariance on the double distance grid,
        % multiplied by the spatial taper in a way that its Fourier
        % transform can be the convolution of the spectral density with the
        % spectral density of the taper, i.e. the expected periodogram
        [Cyy,tyy]=spatmat(ydim,xdim,th,params,xver,dth,tsto);

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
        [Cyy,tyy]=spatmat(ydim,xdim,th,params,xver,dth,tsto);

        % Exploit the symmetry just a tad, which allows us to work with smaller matrices
        q1=fft2(Cyy);
        q4=q1+[q1(:,1) fliplr(q1(:,2:end))];

        % Here is the blurred covariance on the 'complete' grid
        Hh=fftshift(2*real(q4-repmat(fft(Cyy(:,1)),1,NyNx(2)))...
	            -repmat(2*real(fft(Cyy(1,1:end))),NyNx(1),1)...
	            +Cyy(1,1));
        if nargout>2
            % If you ever wanted tyy/Cyy to come out you'll need to unquarter it
            tyy=[fliplr(tyy(:,2:end)) tyy]; tyy=[flipud(tyy(2:end,:)) ; tyy];
            tyy=[zeros(size(tyy,1)+1,1) [zeros(1,size(tyy,2)) ; tyy]];
            Cyy=[fliplr(Cyy(:,2:end)) Cyy]; Cyy=[flipud(Cyy(2:end,:)) ; Cyy];
            Cyy=[zeros(size(Cyy,1)+1,1) [zeros(1,size(Cyy,2)) ; Cyy]];
        end
      case 'efd'
        % Generates a 2*NyNx-1 double grid from which we subsample; applied for
        % calculating the covariance of the DFT over a diagonal in JMATRIX
        ydim=[-NyNx(1)+1:NyNx(1)-1]';
        xdim=[-NyNx(2)+1:NyNx(2)-1] ;

        % Here is the Matern spatial covariance on the double distance grid,
        % multiplied by the spatial taper in a way that its Fourier
        % transform can be the convolution of the spectral density with the
        % spectral density of the taper, i.e. the expected periodogram
        [Cyy,tyy]=spatmat(ydim,xdim,th,params,xver,dth,tsto);

        % Implementing Arthur's fold for the univariate, 2-dimensional case
        numvars=1;
        cbar=reshape(Cyy,[size(Cyy) numvars numvars]);
        foldarr=zeros([NyNx,numvars,numvars]);
        numdims=2;
        for ind=0:numdims-1
          for jnd=0:numdims-1
            res=cbar(ind*NyNx(1)+1:min([(ind+1)*NyNx(1) NyNx(1)*2-1]),...
                     jnd*NyNx(2)+1:min([(jnd+1)*NyNx(2) NyNx(2)*2-1]));
            pres=[zeros(size(res,1)+ind,jnd) [zeros(ind,size(res,2)); res]];
            foldarr=foldarr+pres;
          end
        end
        % To get Cyy out as the size of the data
        Cyy=foldarr;
        Hh=fft2(foldarr);
    end

    % Normalize and vectorize
    Sbar=Hh(:)*prod(dydx)/(2*pi)^2;

    % Should check positivity if you want the non-derivative outputs
    if any(Sbar<0) & isempty(dth) & ~any(tsto); keyboard; end

    % Check Hermiticity of the results for the straight variance case
    if xver==1 & isempty(dth) & ~any(tsto)
        blurcheck(Sbar,params)
        if strcmp('method','ef')
            hermcheck(tyy)
            hermcheck(Cyy)
        end
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

    % Method of computation; default does not assume symmetry
    if isnumeric(method)
      mm=method;
      method=[];
    end
    defval('method','ef')
    % Type of taper, e.g. boxcar, France, Ukraine, Speckle, % Cosine,
    %                     Rectangular, Structured 
    defval('tsto',1)
    % Include taper in the analysis (i.e., calculate the expected periodogram by
    % correctly accounting for variations in sampling); or if 1, treat with
    % square taper
    defval('dth',0)
    % Field size 
    defval('params',[188 233]+randi(20,[1 2]))

    % Number of iterations
    defval('xver',100)
    % Option to use xver input to specify datum of existing calculations to
    % load for sx2 visual summary, provided as string 
    if isstr(xver)
        datum=xver;
        clear xver;
        defval('xver',100);
        % Additions from EGGERS6 for sx2 visual statistics summary require loading
        % MLEOSL batch calculations
        try
            [th,thhats,p,~,~,~,~,~,~,~,momx]=osload(datum,[],'mleosl');
            warning(sprintf('%s %s\n%s %i %s',...
                            'OSLOAD retrieved existing moment parameter data for',datum,...
                            'make sure that taper type',tsto,'corresponds to these data'))
        catch
            % If premade simulations do not exist, generate a suite using
            % MLEOSL defaults for th and p. These values will be retrieved from
            % OSLOAD. (Recall that MLEOSL('DEMO1') will store simulations for
            % today's date only)
            disp(sprintf('Calculating %i simulations for new datum %s',xver,date))
            mleosl('demo1',xver);
            datum=date;
            [th,thhats,p,~,~,~,~,~,~,~,momx]=osload(datum,[],'mleosl');
        end 
    else 
        % BEGIN Figure for boxcar paper -- this will not include the sx2 plots
        %params=[33 33]*3;
        % Some Matern paramters
        th=[1e2 1.5 2.8e3];
        th=[1e2 0.75 2e3];
        th=[1.15 0.65 1.8e3];
        th=[1e4 1.5 5.0e3];
        % Some number of iterations
        % xver=100;
        % Taper type
        % tsto=1;
        % END Figure for boxcar paper

        % Some combinations - SIMULOSL
        p.NyNx=[params(1) params(2)];
        p.dydx=1e3*[1 1];
        p.blurs=Inf;
    end

    % Compare the average periodogram with the blurred spectral density
    
    % What kind of a test are we running?
    if isstr(tsto)
      p.mask=tsto;
    else
      switch tsto
        case 1
          % Boxcar is 0 or 1 or ones, the unit taper
          p.taper=1;
          p.mask='boxcar';
          % The mask parameter is irrelevant, though it might be the anti-taper
        case 2
          % Here's another one
          p.mask='france';
        case 3
          % Here's another one
          p.mask='ukraine';
        case 4
          % Here is the speckled case; default is 66.7% observed
          p.mask='random';
        case 5
          % Here is the % cosine taper
          % TO BUILD
          p.mask=gettaper();
        case 6
          % Here is the rectangle case, which differs from 1 by padding the brims
          % of the field with zeros
          % TO BUILD
          p.mask=gettaper();
        case 7
          % The structured case; for now, this works a lot like the France and
          % Ukraine examples as it is a static mask from a single data sampling
          % pattern extracted from GEBCO 2024 Grid for direct (single beam and
          % multibeam) measurements where C11 is approx at 17S 35.1W and CMN 
          % 12.7S 31.3W for a 4-times decimated set of coordinates from the 
          % original dataset's 15 arc-sec resolution
          % Recall that the amount of data covered by this mask is 62.4%
          % We can take the mask as an arbitrary mask
          geb=load('MYMASKS/gebcostruct.mat');
          % Or we could do this more elaborately with the physical units of the data
          geb=load('MYMASKS/GEBCO2024Grid_maskregion.mat');
          dy=unique(geb.lat_decm-circshift(geb.lat_decm,1));
          dx=unique(geb.lon_decm-circshift(geb.lon_decm,1));
          % Ignore the first output as this is a consequence of using circshift
          if numel(dy)>2; diferm(dy(2),dy(3)); end
          if numel(dx)>2; diferm(dx(2),dx(3)); end
          % Spacing in meters
          p.dydx=[deg2km(dy(2))*1000 deg2km(dx(2))*1000];
          Ny=size(geb.lat_decm,1);
          Nx=size(geb.lon_decm,1);
          p.mask='structured';
          p.taper=logical(geb.I_decm);
          if all(p.NyNx<= [Ny Nx])
            p.taper=p.taper(1:p.NyNx(1),1:p.NyNx(2));
          else
            p.NyNx=[Ny Nx];
          end
        case 8
          % Smoothed France - France with smooth taper applied to exterior boundary
          p.taper=maskit('demo5');
        case 9
          % A random mask
          XY=getMask;
          p.taper=[];
      end
    end
    % Option to shift the mask's center relative to the grid's center?
    % p.mymx=[];
    if isstr(tsto)
         scl=0.8;
         % rotate and shift
         p.mr=180;
         p.mymx=[-10 -8];
        % Generate "mask" only, never mind what the data will be
        [~,I]=maskit(rand(p.NyNx),p,scl);
        % Keep the mask as a taper or the taper as mask, for illustration only
        p.taper=I;
    else
      switch tsto
        case {2,3}
          % To coordinate with the GEBCO "structure" mask,
          scl=1.197;
          p.NyNx=[117 124];
          % AntiFrance: scl=0.585;

          % Generate "mask" only, never mind what the data will be
          [~,I]=maskit(rand(p.NyNx),p,scl);
          % Keep the mask as a taper or the taper as mask, for illustration only
          p.taper=I;
        case 4
          % Generate "(anti-)mask" only, never mind what the data will be
          scl=1-0.667;
          % unless we are here from IAMGwork
          [~,funnam]=star69;
          if strcmp(funnam,'IAMGwork')
            scl=0.7031;
          end
          [~,I]=muckit(rand(p.NyNx),p,scl);
          % Keep the speckles as a taper or the taper as speckles, for
          % illustration only
          p.taper=~I;
          if strcmp(funnam,'IAMGwork')
            p.taper=I;
          end
      end
    end
    % To pass the rigorous test for the residuals, apply a smoothing taper to
    % the outer edges of the field; assuming a rectangle for now, where this
    % works wonders
    if ismember(mm,[1,2,3])
      if tsto==1
        Tx=gettaper(p,'cosine',0.1);
      elseif tsto==4
        % random doesn't suffer edge effects the same way as the contigous
        % windows
        Tx=1;
      elseif ismember(tsto,[7])
        % structured just doesn't look good using the smoothing algorithm and 
        % it does just about as well on edge effects as random; we will smooth
        % the outermost edges anyway
        tp=[];tp.NyNx=p.NyNx;tp.dydx=p.dydx;tp.taper=1;
        Tx=gettaper(tp,'cosine',0.1);
      else
        % we are in a contiguous or structured mask scenario, place smoothing on
        % boundaries
        Tx=gettaper(p,'cosine',[0.9 0.2]);
      end
      p.taper=p.taper.*Tx;
    end


    % Calculate expected periodogram, i.e. the appropriately blurred likelihood
    % Use the exact method via BLUROSY, force p.blurs=-1
    p.blurs=-1;
    % If we want to look at the effect of NOT accounting for sampling
    % irregularities in creating the periodogram, set the taper to 0 following
    % data creation in the analysis step
    if dth==1
        pb=p; pb.taper=1;
    else
        pb=p;
    end
    % Just do the real xver=1 explicitly here
    Sbar=blurosy(th,pb,1,method);
    % Then for what comes next, to simulate data using SGP, force p.blurs=Inf
    p.blurs=Inf;

    % One random one from the sequence will be shown
    randix=randi(xver);
    
    % We be collecting the average periodogram to compare with the expected periodogram
    Sbb=0;
    % The third input in the demo was the number of iterations, fake-called 'xver'
    for index=1:xver
        % Simulate sequentially, collect the expected periodogram, and
        % make the average periodogram
        [Hx,th0,p,k,Hk,Sb,Lb,gane,miy]=simulosl(th,p);
        % Ignore this assessment if we are purposefully messing things up
        if dth~=1
          % Are none of them really bad?
          if max(Sb./Sbar)>20;keyboard; warning('ratio seems bad') ; end
        end
        if index==randix;
            % Keep the special ones for plotting later
            Hxx=Hx; Sbx=Sb;
            % This is the special one that we shall plot later on
            Xk=Sb./Sbar; varibal='X';
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
    % Number of standard deviations for the axis limits and color bars
    nsig=2;

    % First figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    clf
    set(gcf,'InnerPosition',[1 1 450*2 450*1.5])
    [ah,ha,H]=krijetem(subnum(2,3));
    axes(ah(1))
    % imagefnan([1 1],p.NyNx([2 1]),v2s(Hxx,p),[],th(1)*[-3 3]); axis image ij
    h(1)=imagefnan([1 1],p.NyNx([2 1]).*p.dydx([2 1])./1e3,v2s(Hxx,p),[],[min(Hxx) max(Hxx)],[],1);
    axis image 
    t(1)=title(sprintf('field # %i',randix));
    yl(1)=ylabel('northing [km]');
    xl(1)=xlabel('easting [km]');
    %
    [cb(1),xcb(1)]=addcb([],[min(Hxx) max(Hxx)],[min(Hxx) max(Hxx)],'kelicol',[],1);
    cb(1).XTick=linspace(min(Hxx),max(Hxx),5);
    cb(1).XTickLabel=round(cb(1).XTick,0);
    cb(1).Position(1)=0.13;
    cb(1).Position(3)=0.213;
    cb(1).Position(2)=0.55;
    xcb(1).String='$\mathcal{H}(\mathbf{x})$ [km]';
    xcb(1).Interpreter='latex';
    xcb(1).FontSize=10;
    
    axes(ah(2))
    prdgm=log10(v2s(Sbx,p));
    % TODO shown in spectral units
    imagesc(prdgm); 
    axis image xy
    colormap(ah(2),parula)
    t(2)=title(sprintf('periodogram # %i',randix));
    xl(2)=xlabel('wavenumber');
    %
    [cb(2),xcb(2)]=addcb([],[min(prdgm(:)) max(prdgm(:))],...
                            [min(prdgm(:)) max(prdgm(:))],'parula',[],0);
    cb(2).XTick=linspace(min(prdgm(:)),max(prdgm(:)),5);
    cb(2).XTickLabel=round(cb(2).XTick,2);
    cb(2).Position(1)=0.41;
    cb(2).Position(3)=0.213;
    cb(2).Position(2)=0.55;
    xcb(2).String='log power spectral density [km^3]';
    xcb(2).FontSize=8.8;

    % The spatial domain taper
    axes(ah(4))
    if length(p.taper)==1; p.taper=ones(p.NyNx); end
    cm=gray(2);
    h(4)=imagesc([1 p.NyNx(2).*p.dydx(2)./1e3],[1 p.NyNx(1).*p.dydx(1)./1e3],double(p.taper));
    axis image xy
    colormap(ah(4),cm)
    %imagesc(p.taper); axis image
    %colormap(ah(4),gray)
    t(4)=title('taper');
    xl(4)=xlabel('easting [km]');
    yl(4)=ylabel('northing [km]');
    %
    [cb(4),xcb(4)]=addcb([],[0 1],[0 1],cm,[],0);
    cb(4).XTick=[0 1];
    cb(4).XTickLabel=round(cb(4).XTick,2);
    cb(4).Position(1)=0.13;
    cb(4).Position(3)=0.213;
    cb(4).Position(2)=0.06;
    xcb(4).String='$w(\mathbf{x})$';
    xcb(4).Interpreter='latex';
    xcb(4).FontSize=10;

    % The expectation of the periodogram
    axes(ah(3))
    specdens=log10(v2s(Sbar,p));
    % TODO shown in spectral units
    imagesc(specdens); axis image xy
    colormap(ah(3),parula)
    %t(3)=title(sprintf('expectation | %s [%i %0.2g %i]',...
    t(3)=title(sprintf('expectation | %s [%0.2g %g %0.2g]',...
                       '\theta =',th));
    % t(3)=title(sprintf('expectation | %s [%g %g %gx]',...
    %'\theta =',th./[1 1 sqrt(prod(p.dydx))]));
    yl(3)=ylabel('wavenumber'); %moveh(yl(3),435);
    xl(3)=xlabel('wavenumber'); 
    %
    [cb(3),xcb(3)]=addcb([],[min(specdens(:)) max(specdens(:))],...
                            [min(specdens(:)) max(specdens(:))],'parula',[],0);
    cb(3).XTick=linspace(min(specdens(:)),max(specdens(:)),5);
    cb(3).XTickLabel=round(cb(3).XTick,2);
    cb(3).Position(1)=0.691;
    cb(3).Position(3)=0.213;
    cb(3).Position(2)=0.55;
    xcb(3).String='log power spectral density [km^3]';
    xcb(3).FontSize=8.8;
    
    % Then compare with the thing coming out of BLUROSY
    axes(ah(6))
    avgprdm=log10(v2s(Sbb,p));
    imagesc(avgprdm); axis image xy
    t(6)=title(sprintf('average | %i realizations',xver));
    xl(6)=xlabel('wavenumber');
    yl(6)=ylabel('wavenumber');% moveh(yl(6),435);
    [cb(6),xcb(6)]=addcb([],[min(avgprdm(:)) max(avgprdm(:))],...
                            [min(avgprdm(:)) max(avgprdm(:))],'parula',[],0);
    cb(6).XTick=linspace(min(avgprdm(:)),max(avgprdm(:)),5);
    cb(6).XTickLabel=round(cb(6).XTick,2);
    cb(6).Position(1)=0.691;
    cb(6).Position(3)=0.213;
    cb(6).Position(2)=0.06;
    xcb(6).String='log power spectral density [km^3]';
    xcb(6).FontSize=8.8;
    
    
    axes(ah(5))
    rto=v2s(Sbb./Sbar,p);
    imagesc(rto); axis image xy
    t(5)=title(sprintf('(aver / expec), m %4.2f, s %4.2f',...
                       m(end),s(end)));
    xl(5)=xlabel('wavenumber');
    try
        set(ah(5),'clim',m(end)+[-1 1]*nsig*s(end))
    end
    [cb(5),xcb(5)]=addcb([],[min(rto(:)) max(rto(:))],...
                            [min(rto(:)) max(rto(:))],'parula',[],0);
    cb(5).XTick=linspace(min(rto(:)),max(rto(:)),5);
    cb(5).XTickLabel=round(cb(5).XTick,2);
    cb(5).Position(1)=0.41;
    cb(5).Position(3)=0.213;
    cb(5).Position(2)=0.06;
    xcb(5).String='aver / expec';
    xcb(5).FontSize=8.8;
    
    % Clean that sh*t up
    [k,dci,dcn,kx,ky]=knums(p);
    lx=p.NyNx(2);
    ly=p.NyNx(1);
    xlis=unique([1 dci(2) p.NyNx(2)]);
    ylis=unique([1 dci(1) p.NyNx(1)]);
    set([ah(1) ah(4)],'xtick',xlis,...
                      'ytick',ylis,...
                      'xticklabel',[0 ceil(lx/2) lx],...
                      'yticklabel',[0 ceil(ly/2) ly])
                     % 'xticklabel',[-ceil(lx/2) 0 floor(lx/2)],...
                     % 'yticklabel',[-ceil(ly/2) 0 floor(ly/2)])
    set([ah(2) ah(3) ah(5) ah(6)],...
                      'xtick',xlis,...
                      'ytick',ylis,...
                      'xticklabel',{'-\pi','0','\pi'},...
                      'yticklabel',{'-\pi','0','\pi'})
                     % 'xticklabel',[-ceil(lx/2) 0 floor(lx/2)],...
                     % 'yticklabel',[-ceil(ly/2) 0 floor(ly/2)])
    set(ah,'FontSize',9)
    longticks(ah)
    set(ah([3 6]),'YAxisLocation','right')
    set(t,'FontSize',10,'FontWeight','normal')
    movev(cb(1:3),-0.03)
    keyboard
    % serre(H,1,'across')
    %serre(H',0.5,'down')
    %serre([H' cb'],0.75,'down')
    % movev(t,-p.NyNx(1)/15)
    [bh,txh]=label(ah,'ll',9,0);
    if isstr(tsto)
      fname=sprintf('demo_2_%3.3i_%3.3i_%s-%3.3i_fld',p.NyNx,tsto,xver);
    else
      fname=sprintf('demo_2_%3.3i_%3.3i_%3.3i-%3.3i_fld',p.NyNx,tsto,xver);
    end
    if dth==1
        fname=append(fname,'cc1');
    end
    set(gcf,'Renderer','painters')
    figdisp([],fname,[],1)

    % Second figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    clf
    set(gcf,'InnerPosition',[1 1 450*2 450*1.5])
    [ah,ha,H]=krijetem(subnum(2,3));
    df=2;
    % Craft some labels
    xll=[0 3*2*df];
    xlls=[xll(1):df:xll(2)];
    xstr=sprintf('quadratic residual %s',varibal);
    xstr2=sprintf('quadratic residual 2%s',varibal);
    cax=[0 3*df];
    cborien='vert';

    axes(ah(1))
    % Literally follow MLECHIPLOS, to a point, spin off function later on
    [bdens,c]=hist(2*Xk,5*round(log(prod(p.NyNx))+1));
    % Plot the histogram as a bar graph
    bdens=bdens/indeks(diff(c),1)/prod(p.NyNx);
    bb=bar(c,bdens,1);
    % Each of the below should be df/2
    t(1)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
	 varibal,nanmean(Xk),...
	 varibal,nanvar(Xk)));
    set(bb,'FaceC',grey)
    hold on
    % Plot the ideal chi-squared distribution
    refs=linspace(0,max(2*Xk),100);
    hold on
    plot(refs,chi2pdf(refs,df),'Linew',0.5,'Color','k')
    hold off
    longticks(ah(1))
    maxi=0.25;
    ylls=[0:0.1:maxi*(4/df)];
    ylim([0 maxi*(4/df)])
    
    xlim(xll)
    set(ah(1),'xtick',xlls)
    xl(1)=xlabel(xstr2); 
    yl(1)=ylabel('probability density');
    axis square
    movev(t(1),maxi*(4/df)/20)

    axes(ah(2))
    % This should be a straight line folks since the ratio is 1/2 chi^2_2
    h=qqplot(2*Xk,makedist('gamma','a',df/2,'b',2)); axis equal;

    axis image; box on
    longticks(ah(2))
    set(h(1),'MarkerEdge','k')  
    set(h(3),'LineStyle','-','Color',grey)
    % Plot the one-to-one line
    hold on
    xh=get(h,'Xdata'); xh=xh{1};
    yh=get(h,'ydata'); yh=yh{1};
    qq0=plot([0 xll(2)],[0 xll(2)],'k');
    delete(h(1))
    % Replot for more control
    plot(xh,yh,'LineStyle','none','Marker','+','MarkerFaceColor','k',...
         'MarkerEdgeColor','k','MarkerSize',2);
    top(h(3),ah(2))
    delete(get(ah(2),'ylabel'));
    delete(get(ah(2),'title'));
    xlim(xll); ylim(xll)
    xl(2)=xlabel(sprintf('predicted 2%s',varibal));
    yl(2)=ylabel(sprintf('observed 2%s',varibal));
    set(gca,'xtick',xlls,'ytick',xlls)

    allg=~isinf(Xk); sum(allg);

    % Test for departure of chi-squaredness via the "magic" parameter which
    % This is the same as what comes out of LOGLIOSL etc
    magx=nanmean([Xk-df/2].^2);
    neem='\xi'; neem='s_X^2';
    % Do the test whether you accept this as a good fit, or not
    vr=8/length(k(~~k));
    [a,b,c]=normtest(magx,1,vr,0.05);
    % disp(sprintf('NORMTEST %i %5.3f %i',a,b,round(c)))
    if a==0; stp='accept'; else stp='reject'; end
    t(2)=title(sprintf('%s = %5.3f  8/K = %5.4f  %s  p = %5.2f',...
                       neem,magx,vr,stp,b));
    movev(t(2),xll(2)/40)
    
    axes(ah(3))
    imagefnan([xlis(1) ylis(end)],[xlis(end) ylis(1)],...
              v2s(Xk,p),'gray',cax,[],1); axis image square
    set(ah(3),'xtick',xlis,'XtickLabel',[-lx/2 0 lx/2],...
        'ytick',ylis,'YtickLabel',[-ly/2 0 ly/2]);
    longticks(ah(3))
    xl(3)=xlabel('wavenumber index'); 
    yl(3)=ylabel('wavenumber index');
    t(3)=title(sprintf('%s | %i x %i | %s = [%0.2g %g %0.2g]',p.mask,p.NyNx,...
                       '\theta',th));
    % t(3)=title(sprintf('%s | %i x %i | %s = [%g %g %gx]',p.mask,p.NyNx,n'\theta',...
                       %th./[1 1 sqrt(prod(p.dydx))]));
    movev(t(3),ylis(end)/20)
    
    % Check whether we chose to include the sx2 visual summaries from a MLEOSL
    % batch run, and if so, make additional figures
    if ~exist('datum','var')
        delete(ah(4:6))

        % Cosmetics from EGGERS6
        axes(ah(3))
        [cb(1),xcb(1)]=addcb(cborien,cax,cax,'gray',df,1);
        axes(cb(1))
        set(xcb(1),'string',xstr)
        set(cb(1),'YAxisLocation','right')
        set(cb(1),'position',...
               [getpos(ah(3),1)+getpos(ah(3),3)*1.1 getpos(ah(3),2) getpos(cb(1),3) getpos(ah(3),4)])
        shrink(cb(1),1,1.15)

        moveh(ah(1),-0.03); moveh(ah(3),0.03)
        moveh(ah(1:3),-0.03)
        movev([ah(1:3) cb],-0.2)

        % Control all axes and font sizes at the same time
        set([ah(1:3) cb],'FontSize',8)
        set(t,'FontSize',9,'FontWeight','normal')
    else
        axes(ah(4))
        mx3=momx(:,3);
        % Take out the Inf which may occur at zero wavenumber
        [bdens,c]=hist(mx3,4*round(log(length(mx3))));
        bdens=bdens/indeks(diff(c),1)/length(mx3);
        bb=bar(c,bdens,1);
        varibal='s_X^2';
        % Each of the below should be df/2
        t(4)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
                   varibal,nanmean(mx3),...
                   varibal,nanvar(mx3)));
        set(bb,'FaceC',grey)
        hold on
        sfax=4;
        k=knums(p); varpred=8/[length(k(~~k))];
        xll=[1-sfax*sqrt(varpred) 1+sfax*sqrt(varpred)];
        xls=linspace(1-sfax*sqrt(varpred),1+sfax*sqrt(varpred),sfax+1);
        refs=linspace(xll(1),xll(2),100);
        plot(refs,normpdf(refs,1,sqrt(varpred)),'Linew',1,'Color','k')
        axis square
        xlim(xll)
        set(ah(4),'xtick',xls,'xticklabel',round(xls*100)/100)
        longticks(ah(4))
        xl(1)=xlabel(varibal);
        yl(1)=ylabel('probability density');
        %ylim([0 11])
        hold off

        % Save this position as it gets messed up afterwards
        posh=getpos(ah(4));
        
        axes(ah(5))
        h=qqplot(mx3,makedist('normal','mu',1,'sigma',sqrt(varpred)));
        axis equal; box on
        set(h(1),'MarkerE','k')
        set(h(3),'LineS','-','Color',grey)
        % Extend the line to the full axis
        xll=[1-sfax*sqrt(varpred) 1+sfax*sqrt(varpred)];
        
        hold on
        xh=get(h(3),'xdata');
        yh=get(h(3),'ydata');
        % On the right
        h(4)=plot([xh(2) xll(2)],...
              [yh(2) yh(2)+[yh(2)-yh(1)]/[xh(2)-xh(1)]*[xll(2)-xh(2)]]);
        % On the left
        h(5)=plot([xll(1) xh(1)],...
              [yh(1)+[yh(2)-yh(1)]/[xh(2)-xh(1)]*[yh(2)-xll(2)] yh(1)]);
        set(h(4:5),'LineS','-','Color',grey)
        hold off
        top(h(3),ah(2))
        delete(get(ah(5),'ylabel'));
        delete(get(ah(5),'title'));
        sfax=3;
        xls=linspace(1-sfax*sqrt(varpred),1+sfax*sqrt(varpred),(2*sfax+1));
        xlc=num2cell(round(xls*10)/10);
        xlim(xll); ylim(xll)
        xlc{1}='';
        xlc{3}='';
        xlc{5}='';
        xlc{7}='';
        set(ah(5),'xtick',xls,'xticklabel',xlc,...
                  'ytick',xls,'yticklabel',xlc)
        longticks(ah(5))
        xl(2)=xlabel(sprintf('predicted %s',varibal));
        yl(2)=ylabel(sprintf('observed %s',varibal));
        axis square

        axes(ah(6))
        sfax=2;
        if th0(1)>=1000
            cax2=[-sfax sfax]*sqrt(th0(1))/1000;
            imagefnan([xlis(1) xlis(end)],[xlis(end) xlis(1)],...
              reshape(Hx-mean(Hx),p.NyNx)/1000,...
              'gray',cax2,[],1); axis image
            xstr3='simulated field [km]';
        else
            cax2=[-sfax sfax]*sqrt(th0(1));
            imagefnan([xlis(1) xlis(end)],[xlis(end) xlis(1)],...
              reshape(Hx-mean(Hx),p.NyNx),...
              'gray',cax2,[],1); axis image
            xstr3='simulated field [m]';
        end
        cborien='vert';
        set(ah(6),'xtick',xlis,'xticklabel',[-lx/2 0 lx/2],...
              'ytick',ylis,'yticklabel',[-ly/2 0 ly/2])
        xl(6)=xlabel('position index');
        yl(6)=ylabel('position index');
        longticks(ah(6))

        % Cosmetics from EGGERS6
        axes(ah(3))
        [cb(1),xcb(1)]=addcb(cborien,cax,cax,'gray',df,1);
        axes(cb(1))
        set(xcb(1),'string',xstr)
        set(cb(1),'YAxisLocation','right')
        set(cb(1),'position',...
               [getpos(ah(3),1)+getpos(ah(3),3)*1.1 getpos(ah(3),2) getpos(cb(1),3) getpos(ah(3),4)])
        shrink(cb(1),1,1.15)

        axes(ah(6))
        if th0(1)>=1000 
            [cb(2),xcb(2)]=addcb(cborien,cax2,cax2,'gray',sqrt(th0(1))/1000,1);
        else
            [cb(2),xcb(2)]=addcb(cborien,cax2,cax2,'gray',sqrt(th0(1)),1);
        end
        axes(cb(2))
        set(xcb(2),'string',xstr3)
        cblabs=sprintfc('%0.1f',str2num(get(cb(2),'YTickLabel')));
        set(cb(2),'YAxisLocation','right','position',...
           [getpos(ah(6),1)+getpos(ah(6),3)*1.1 getpos(ah(6),2) getpos(cb(2),3) getpos(ah(6),4)],...
           'YTickLabel',cblabs);
        shrink(cb(2),1,1.15)
        moveh(ah([1 4]),-0.03);moveh([ah([3 6]) cb],0.03)
        moveh([ah cb],-0.03)

        % Control all axes and font sizes at the same time
        set([ah cb],'FontSize',8)
        set(t,'FontSize',7,'FontWeight','normal')
        pah3b=getpos(ah(3));
        serre(H',0.25,'down')
        pah3a=getpos(ah(3));
        movev(cb(1),pah3a(2)-pah3b(2));
    end

    if isstr(tsto)
      fname=sprintf('demo_2_%3.3i_%3.3i_%s-%3.3i_chi',p.NyNx,tsto,xver);
    else
      fname=sprintf('demo_2_%3.3i_%3.3i_%3.3i-%3.3i_chi',p.NyNx,tsto,xver);
    end
    if dth==1
        fname=append(fname,'cc1');
    end
    figdisp([],fname,[],1)
    %figdisp([],sprintf('demo_2_%3.3i_%3.3i_%3.3i-%3.3i_chi',p.NyNx,tsto,xver),[],1)

    % Third figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3)
    clf
    H=errorbar(1:xver,m,-nsig*s,nsig*s);
    try
        set(getkids(H,1),'LineWidth',0.5)
        set(getkids(H,2),'Color',grey); 
    catch
        H.LineWidth=0.5;
        H.Color=grey;
    end
    longticks(gca,2)
    axis tight; grid on; xlabel('sample size'); 
    ylabel(sprintf('average / expected periodogram | %i sigma',nsig))
    xlim([0 xver+1])
    set(gca,'xtick',[1 10:10:xver])
    shrink(gca,1,1.5)
    ylim([min(m-s*2) max(m+2*s)]*1.1)

    hold on
    plot(1:xver,m,'LineWidth',1,'Color','b')
    plot(1:xver,1.96./sqrt(1:xver)+m,'LineWidth',1,'Color','k')
    hold off
    %tt=title(sprintf('%s | %i x %i | %s = [%g %g %gx]',p.mask,p.NyNx,'\theta',...
    %                 th./[1 1 sqrt(prod(p.dydx))]));
    tt=title(sprintf('%s | %i x %i | %s = [%0.2g %g %0.2g]',p.mask,p.NyNx,...
                     '\theta',th));
    movev(tt,range(ylim)/30)
    if isstr(tsto)
      fname=sprintf('demo_2_%3.3i_%3.3i_%s-%3.3i_std',p.NyNx,tsto,xver);
    else
      fname=sprintf('demo_2_%3.3i_%3.3i_%3.3i-%3.3i_std',p.NyNx,tsto,xver);
    end
    if dth==1
        fname=append(fname,'cc1');
    end
    figdisp([],fname,[],1)
    %figdisp([],sprintf('demo_2_%3.3i_%3.3i_%3.3i-%3.3i_std',p.NyNx,tsto,xver),[],1)
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
    tolex=-(log10(max(Sbar1))-11);
    diferm(Sbar1,Sbar2,tolex);
    diferm(tyy1,tyy2); diferm(Cyy1,Cyy2);
    diferm(Sbar1,Sbar3,tolex);
    diferm(tyy1,tyy3); diferm(Cyy1,Cyy3);
    diferm(Sbar1,Sbar4,tolex);
    diferm(tyy1,tyy4); diferm(Cyy1,Cyy4);
    diferm(Sbar2,Sbar3,tolex);
    diferm(tyy2,tyy4); diferm(Cyy2,Cyy3);
    diferm(Sbar2,Sbar4,tolex);
    diferm(tyy2,tyy4); diferm(Cyy2,Cyy4);
    diferm(Sbar3,Sbar4,tolex);
    diferm(tyy3,tyy4); diferm(Cyy3,Cyy4);
    % subplot(221); imagesc(log10(v2s(Sbar1,p))); axis square
elseif strcmp(th,'demo4')
    % Look at a row of the expected and averaged periodograms of the interior
    % and exterior of an irregular boundary on a small grid.
    th=[1 0.5 2];
    p=[];p.NyNx=[9 9];p.dydx=[1 1];p.blurs=-1;p.mask='newjersey';
    [~,I]=maskit(ones(p.NyNx),p,[]);
    p.taper=I;
    Sbar=blurosy(th,p,1,'ef');
    p2=p;p2.taper=~I;
    Sbar2=blurosy(th,p2,1,'ef');

    % collect the average periodograms
    Sbb=0;
    Sbb2=0;
    % The third input in the demo was the number of iterations, fake-called 'xver'
    for index=1:xver
        % Simulate sequentially, collect the expected periodogram, and
        % make the average periodogram
        p.blurs=Inf;
        [~,~,~,~,~,Sb]=simulosl(th,p);
        [~,~,~,~,~,Sb2]=simulosl(th,p2);
        % Are none of them really bad?
        if max(Sb./Sbar)>20;keyboard; warning('ratio seems bad') ; end
        if max(Sb2./Sbar2)>20;keyboard; warning('ratio seems bad') ; end
        % Collect the average, watch the growth later
        Sbb=Sbb+Sb;
        Sbb2=Sbb2+Sb2;
        % Mean of the evolving ratio over the wavenumbers
        m(index)=mean(Sbb/index./Sbar);
        % Mean of the evolving standard deviation of the mean ratio over the wavenumbers
        s(index)=std(Sbb/index./Sbar);
    end 
    % Mean ratio over the realizations
    Sbb=Sbb/xver;
    Sbb2=Sbb2/xver;

    figure();
    % this was nice
    %set(gcf,'InnerPosition',[1 1 1145 273])
    % but we will use this to accomodate the colorbar
    set(gcf,'InnerPosition',[1 1 1145 300])
    [ah,ha]=krijetem(subnum(1,4));
keyboard
    % expected, interior
    axes(ah(1))
    imagesc(log10(v2s(Sbar,p))); axis image                                     
    hold on
    yl(1)=ylabel('$k_y$','Interpreter','latex');%'wavenumber index');
    xl(1)=xlabel('$k_x$','Interpreter','latex');%'wavenumber index');  
    cb(1)=colorbar; cb(1).Label.String='log, [units^2 length]';
    cb(1).Location='south';
    cb(1).Ticks=linspace(min(caxis),max(caxis),5);
    cb(1).TickLabels=round(cb(1).Ticks,1);
    cb(1).Position(2)=0.03;
    cb(1).Position(4)=0.03;
    cb(1).Position(3)=0.1436;
    cb(1).AxisLocation='out';
    cb(1).TickDirection='out';
    t(1)=title('$\bar{\mathcal{S}}_{\boldmath{\theta}}(\mathbf{k})$','Interpreter','latex');
    inax(1)=axes('position',[0.234 0.238 0.05 0.20]);
    imagesc(I)
    colormap(inax(1),gray)

    % averaged, interior
    axes(ah(2))
    imagesc(log10(v2s(Sbb,p))); axis image                                     
    t(2)=title('$\overline{|H(\mathbf{k})|^2}$','Interpreter','latex');
    yl(2)=ylabel('$k_y$','Interpreter','latex');%'wavenumber index');
    xl(2)=xlabel('$k_x$','Interpreter','latex');%'wavenumber index');                                           
    cb(2)=colorbar; cb(2).Label.String='log, [units^2 length]';
    cb(2).Location='south';
    cb(2).Ticks=linspace(min(caxis),max(caxis),5);
    cb(2).TickLabels=round(cb(2).Ticks,1);
    cb(2).Position(2)=0.03;
    cb(2).Position(4)=0.03;
    cb(2).Position(3)=0.1436;
    cb(2).AxisLocation='out';
    cb(2).TickDirection='out';
    
    % expected, exterior
    axes(ah(3))
    imagesc(log10(v2s(Sbar2,p2))); axis image                                     
    hold on
    t(3)=title('$\bar{\mathcal{S}}_{\boldmath{\theta}}(\mathbf{k})$','Interpreter','latex');
    xl(3)=xlabel('$k_x$','Interpreter','latex');%'wavenumber index'); 
    yl(3)=ylabel('$k_y$','Interpreter','latex');%'wavenumber index');
    cb(3)=colorbar; cb(3).Label.String='log, [units^2 length]';
    cb(3).Location='south';
    cb(3).Ticks=linspace(min(caxis),max(caxis),5);
    cb(3).TickLabels=round(cb(3).Ticks,1);
    cb(3).Position(2)=0.03;
    cb(3).Position(4)=0.03;
    cb(3).Position(3)=0.1436;
    cb(3).AxisLocation='out';
    cb(3).TickDirection='out';
    inax(2)=axes('position',[0.646 0.238 0.05 0.20]);
    imagesc(~I)
    colormap(inax(2),gray)

    % averaged, exterior
    axes(ah(4))
    imagesc(log10(v2s(Sbb2,p2))); axis image                                     
    t(4)=title('$\overline{|H(\mathbf{k})|^2}$','Interpreter','latex');
    xl(4)=xlabel('$k_x$','Interpreter','latex');%'wavenumber index');                                           
    yl(4)=ylabel('$k_y$','Interpreter','latex');%'wavenumber index');
    cb(4)=colorbar; cb(4).Label.String='log, [units^2 length]';
    cb(4).Location='south';
    cb(4).Ticks=linspace(min(caxis),max(caxis),5);
    cb(4).TickLabels=round(cb(4).Ticks,1);
    cb(4).Position(2)=0.03;
    cb(4).Position(4)=0.03;
    cb(4).Position(3)=0.1436;
    cb(4).AxisLocation='out';
    cb(4).TickDirection='out';

    % clean up                                                        
    [k,dci,dcn,kx,ky]=knums(p);                                                 
    lx=p.NyNx(2);                                                               
    ly=p.NyNx(1);                                                               
    %xlis=unique([1 dci(2) p.NyNx(2)]);                                          
    %ylis=unique([1 dci(1) p.NyNx(1)]);                                          
    xlis=unique([1 p.NyNx(2)]);                                          
    ylis=unique([1 p.NyNx(1)]);                                          
    set([ah(1) ah(2) ah(3) ah(4)],...                                           
                      'xtick',xlis)                                          
                      %'xticklabel',{'-\pi','0','\pi'})
    set([ah(1) ah(2) ah(3) ah(4)],...                                           
                      'ytick',ylis)                                         
                      %'yticklabel',{'-\pi','0','\pi'})                          
    %set([ah(2) ah(3)],'yticklabel',{'',''})
    set([inax(1) inax(2)],'ytick',{},'xtick',{})
    %set(ah(4),'YAxisLocation','right')                                      
    set(ah,'FontSize',9)                                                        
    longticks(ah)                                                               
    set(t,'FontSize',10,'FontWeight','normal')                                  
    movev(t,-p.NyNx(1)/50)
    movev([ah cb inax],0.07)
    movev(cb,0.03)
    movev(inax,0.001)
    [bh,txh]=label(ah,'ul',9,0);
    uistack(inax,'top')
    set(gcf,'Renderer','painters')
    keyboard

    % save
    fname=sprintf('demo_4_%3.3i_%3.3i_%s',p.NyNx,p.mask);
    figdisp([],fname,[],1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cyy,t]=spatmat(ydim,xdim,th,params,xver,dth,tsto)
% [Cyy,t]=spatmat(ydim,xdim,th,params,xver,dth,tsto)
%
% Returns the modified spatial covariance whose Fourier transform is the blurred
% spectrum after spatial data tapering, i.e. the expected periodogram. No use to
% return the autocovariance of the applied spatial taper which is an essential
% part of this operation...  never need it explicitly, used in SIMULOSL and
% LOGLIOSL via the intermediary of BLUROSY.

% Dimensions of the original grid
NyNx=params.NyNx;
dydx=params.dydx;
if ~isfield(params,'taper')
    params.taper=[];
end

% Specify the spatial taper EXPLICITLY
if prod(size(params.taper))>1
    % Now compare with the other mechanisms
    % Completely general windows where 1 means you are taking a sample
    % See Arthur's note for more general windows, use IFF2/FFT2 you need, see
    % ~/POSTDOCS/ArthurGuillaumin/CodeArthur/NonParametricEstimation/Periodogram.m
    % $MFILES/retired/QR?.m/map_*.m/whittle_*/expected_* etc
    Tx=double(params.taper);

    % If you are here with 'efs' the taper is explicit AND not
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
        % top way, but leave the explicit way for illustration - this is
        % the long step
        t=xcorr2(Tx);
	% Add a row of zeros here
	t=[zeros(size(t,1)+1,1) [zeros(1,size(t,2)) ; t]];
    end
    % Now normalize the cross-correlations at the end
    t=t/sum(sum(Tx.^2));

    % disp(sprintf('\nStill can optimize\n'))
    % Here too should use FFT where we can, see COMPUTE_KERNELS
    % internally and below, I would imagine that's just the same thing
    % inside Arthur's PERIODOGRAM object
    % kernel1 = ifft2(abs(fft2(Tx,2*p.NyNx(1),2*NyNx(2))).^2/prod(p.NyNx));
    % kernel1 = kernel1(1:p.NyNx(1),1:p.NyNx(2));
    % kernel2 = ifft2(abs(fft2(fliplr(Tx),2*p.NyNx(1),2*p.NyNx(2))).^2/prod(p.NyNx));
    % kernel2 = kernel2(1:p.NyNx(1),1:p.NyNx(2));
    % Then feed that into COMPUTE_EXPECTATION_FROM_KERNELS

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
        if ~any(tsto)
          % Fix the rim by adding zeroes top and left
          t2=[zeros(size(t2,1)+1,1) [zeros(1,size(t2,2)) ; t2]];
        end
        % Check the difference between these two implementations,
        % all checked for even/odd/method combinations on 2/24/2023
        if all(size(t)==NyNx)
            % Then the doubling inside this block needs to be undone
            t2=t2(NyNx(1)+1:end,NyNx(2)+1:end);
        end
        diferm(t,t2,9);
    end
end

% Use the tsto input argument for providing frequency index offsets required
% for calculating the JMATRIX per-diagonal method (3; Eq. 17 of VoWE)
if ~isempty(tsto) & size(tsto,2)==2
    m1=tsto(1); m2=tsto(2);
    normalization_factor=prod(NyNx);
    two_n=NyNx.*2-1;
    a=exp(2i*pi*m1/NyNx(1).*(NyNx(1):-1:1))';
    b=reshape(exp(2i*pi*m2/NyNx(2).*(NyNx(2):-1:1))',1,NyNx(2));
    c=a*b;
    if isfield(params,'taper') & numel(params.taper)>1
        g=params.taper;
        sg=sum(g.^2,"all");
        g2=g.*c; 
    else
        g=ones(NyNx);
        sg=numel(g);
        g2=c;
    end
    % We can speed up the calculation sometimes using gpuArray calculations
    % below, however this only pays off for grids larger than about 2500 
    % elements when not using parfor and larger than 6000 elements when using 
    % parfor, which we should
    mkgpuarr=1 & normalization_factor>6e3;
    if mkgpuarr
        g=gpuArray(g);
        g2=gpuArray(g2);
    end
    % The next two lines are the most expensive part of calculating JMATRIX
    % method 3, so we provide gpuArrays when it will help
    f =fftn(g,two_n).*conj(fftn(g2,two_n));
    cg=ifftn(f);
    if mkgpuarr
        cg=gather(cg);
    end
    t =cg./sg;
end

% Here is the distance grid, whose size depends on the input
y=sqrt(bsxfun(@plus,[ydim*dydx(1)].^2,[xdim*dydx(2)].^2));

% The modified spatial covariance
Cyy=maternosy(y,th,dth);
if ~isempty(tsto) & size(tsto,2)==2
    Cyy=ifftshift(Cyy);
end
Cyy=Cyy.*t;
