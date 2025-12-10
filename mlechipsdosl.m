function varargout=mlechipsdosl(Hk,thhat,scl,params,stit,ah,unts)
% [a,mag,ah,ah2,cb,ch,spt,cll]=MLECHIPSDOSL(Hk,thhat,scl,params,stit,ah,unts)
%
% Makes a FOUR-panel plot of the quadratic residuals and their
% interpretation for a Matern likelihood model of a single-field
% Fourier-domain data patch, as well as plots of the actual power spectral
% density of that data patch and its predicted power spectral density based
% on the model. These residuals should have a chi-squared distribution, so
% the produced figure is useful in determining how well the Matern
% likelihood model fits the data.
% 
% The figure produced consists of four panels:  
%    (TL) Histogram and qq-plot of the chi-squared residuals
%    (TR) Predicted 2D power spectral density based on "thhat"
%    (BL) 2D chi-squared residuals
%    (BR) "True" power spectral density with contours from (TR) overlain
%
% INPUT:
%
% Hk         The Fourier-domain data, e.g. from SIMULOSL, or, if all reals, we
%            will take it that the data provided is in the space-domain and we
%            will be performing differencing of the data for comparison of the
%            residuals (e.g., Biometrika 2019)
% thhat      The evaluated scaled Matern parameter vector
% scl        The scale factors
% params     The structure with the fixed parameters from the experiment
% stit       Title for the overall figure [defaulted]
% ah         A quartet of axis handles [defaulted]
% unts       A unit string for the space variable    
%
% OUTPUT:
%
% a           1 or 0 whether the model was accepted or not
% mag         The "magic" parameter which I've come to call s_X^2
% ah,ah1      Main and secondary axis (panels 1-4) handles
% cb          Colorbar (panel 3) handle
% ch          Contour (panels 2 and 4) handles
% spt         Supertitle handle
% cll         Handles to contour labels
%
% SEE ALSO:
%
% MLECHIPLOS, MLEPLOS, EGGERS8
%
% NOTE: 
%
% Maybe should integrate MLECHIPLOS into this one.
%
% EXAMPLES:
%
% p=[]; p.NyNx=[64 64]; p.dydx=[1e4 1e4]; p.kiso=NaN;
% th=[1e6 2.5 2e4];
%
% % Try with a smoothing taper:
% Tx=gettaper(p,'cosine',0.1);
%
% % Sometimes works well:
% p.blurs=Inf; p.taper=1;
% [Hx1,th,p,~,Hk1]=simulosl(th,p); 
% p.blurs=-1; p.taper=Tx;
% [thhat1,~,~,scl1,~,~,Hk1]=mleosl(Hx1.*p.taper(:),[],p);scl1(1)=1;
% mlechipsdosl(Hk1,thhat1,scl1,p)
%
% % Works great:
% p.blurs=-1; p.taper=1;
% [Hx2,~,~,~,Hk2]=simulosl(th,p);
% p.blurs=-1; p.taper=Tx;
% [thhat2,~,~,scl2,~,~,Hk2]=mleosl(Hx2.*p.taper(:),[],p);scl2(1)=1;
% mlechipsdosl(Hk2,thhat2,scl2,p)
%
% % Not as expected:
% th(2)=4.5;
% p.blurs=Inf; p.taper=1;
% [Hx3,~,~,~,Hk3]=simulosl(th,p);
% p.blurs=-1; p.taper=Tx;
% [thhat3,~,~,scl3,~,~,Hk3]=mleosl(Hx3.*p.taper(:),[],p);scl3(1)=1;
% mlechipsdosl(Hk3,thhat3,scl3,p)
%
% % Very much not as expected:
% th(2)=6.5;
% p.blurs=Inf; p.taper=1;
% [Hx4,~,~,~,Hk4]=simulosl(th,p);
% p.blurs=-1; p.taper=Tx;
% [thhat4,~,~,scl4,~,~,Hk4]=mleosl(Hx4.*p.taper(:),[],p);scl4(1)=1;
% mlechipsdosl(Hk4,thhat4,scl4,p)
% 
% mlechipsdosl('demo1') % Quartzite thin section 
% mlechipsdosl('demo2') % Atlantic bathymetry
% mlechipsdosl('demo3') % Fault roughness
% mlechipsdosl('demo4') % Venus topography
% mlechipsdosl('demo5') % Sea surface height anomaly 
%
% Last modified by gleggers-at-princeton.edu, 04/17/2014
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2018
% Last modified by olwalbert-at-princeton.edu, 06/19/2025

% Some defaults
defval('stit','Chi-squared residuals')
% This prevents the passing of any prior existing handles...
defval('ah',krijetem(subnum(2,2)))
defval('unts','km')
       
if ~isstr(Hk)
    if isreal(Hk)
        % If we provided a real data vector, this is spatial data and we want to
        % investigate whether differencing the data improves the residuals. We will
        % calculate differenced data first, then the corresponding autocovariance
        % sequence for the differenced process.
        shat=nanstd(Hk(:)); 
        Hxd=v2s(Hk./shat,params);
        pd=params;
        % Differencing matrices
        Dy=sparse(diag(ones(prod(pd.NyNx(1))-1,1),1)-diag(ones(prod(pd.NyNx(1)),1),0));
        Dx=sparse(diag(ones(prod(pd.NyNx(2))-1,1),1)-diag(ones(prod(pd.NyNx(2)),1),0));
        Hxd=Dy*Hxd;
        % Again, transpose
        Hxd=v2s(Hxd,pd)';
        Hxd=(Dx*Hxd)';
        Hxd=Hxd(2:end-1,2:end-1);Hxd=Hxd(:);
        pd.NyNx=pd.NyNx-[2 2];
        % That should be the same as this (almost...?)
        % Hxd=v2s(Hk./shat,params);
        % DHx=Hxd(2:end,:)-Hxd(1:end-1,:);
        % DDHx=DHx(:,2:end)-DHx(:,1:end-1);
        % Hxd=DDHx;
        % pd=params;pd.NyNx=size(Hxd);Hxd=Hxd(:);
        
        % Now we have the spectral data
        if numel(pd.taper)>1
            Hk=tospec(Hxd,pd)/(2*pi)/sqrt(sum(pd.taper(:).^2))*sqrt(prod(pd.NyNx));
        else
            Hk=tospec(Hxd,pd)/(2*pi);
        end
        
        % Form the autocovariance sequence of the spatially differenced data, much of
        % this is pulled directly from BLUROSY method 'ef' and SPATMAT subfunction
        NyNx=pd.NyNx;
        dydx=pd.dydx;
        ydim=[-NyNx(1):NyNx(1)-1]';
        xdim=[-NyNx(2):NyNx(2)-1] ;
        
        % The autocovariance sequence for the differenced data, a la Biometrika 2019.
        % Note that this equation held for 1-D; is it still appropriate for our
        % isotropic model?
        y0 =sqrt(bsxfun(@plus,[ydim*dydx(1)].^2,[xdim*dydx(2)].^2));
        yp1=sqrt(bsxfun(@plus,[(ydim-1)*dydx(1)].^2,[xdim*dydx(2)].^2));
        ym1=sqrt(bsxfun(@plus,[(ydim+1)*dydx(1)].^2,[xdim*dydx(2)].^2));
        xp1=sqrt(bsxfun(@plus,[ydim*dydx(1)].^2,[(xdim+1)*dydx(2)].^2));
        xm1=sqrt(bsxfun(@plus,[ydim*dydx(1)].^2,[(xdim-1)*dydx(2)].^2));
        
        scl2=scl;scl2(1)=1; ths=thhat.*scl2;
        % No, it is not. This results in a blurred spectrum with an erroneous
        % lemniscate geometry.
        % Cyd=2*maternosy(y0,ths)-...
        %       maternosy(yp1,ths)-maternosy(ym1,ths);
        
        % We should probably adapt the modifications of the lags as a quincunx.
        Cyd=4*maternosy(y0,ths)-...
            maternosy(yp1,ths)-maternosy(ym1,ths)-...
            maternosy(xp1,ths)-maternosy(xm1,ths);
        
        % We can apply the taper in space as usual, taking care to trim the outer-most 
        % edge to match the new data size.
        if numel(pd.taper)>1 
            pd.taper=pd.taper(2:end-1,2:end-1);
            Tx=double(pd.taper); 
        else 
            Tx=ones(NyNx);
        end
        t=xcorr2(Tx);
        t=[zeros(size(t,1)+1,1) [zeros(1,size(t,2)) ; t]];
        t=t/sum(sum(Tx.^2));
        Cyd=Cyd.*t;
        
        % Calculate the blurred spectrum for the differenced process
        Hh=fftshift(realize(fft2(ifftshift(Cyd))));
        Hh=Hh(1+mod(NyNx(1),2):2:end,1+mod(NyNx(2),2):2:end);
        Sb=Hh(:)*prod(dydx)/(2*pi)^2;
        
        % We will need a set of wavenumbers for later and things will be easier if we
        % reassign params.
        params=pd;
        [kk,~,~,kx,ky]=knums(pd); 
        k=kk(:);

        % See Vanmarcke 2010 differencing
    else
        % Same as it always was:
        % The following is as in MLECHIPLOS, need to force SIMULOSL to return theoretical Sb
        % note that in MLECHIPLOS we still had a number of alternatives to compute the same things
        params.blurs=-1;

        % So now you have found a solution and you evaluate it at thhat
        [~,~,~,k,~,Sb,Lb]=simulosl(thhat.*scl,params,0);
        % Get just a little more information
        [kk,~,~,kx,ky]=knums(params); 
    end
    
    % Degrees of freedom for 2X system (assuming a single field)
    df=2;
    % Get the quadratic residuals, any of a number of ways

    % Multiply to obtain a variable which should follow the rules 
    %Zk=[Hk(:)./Lb(:)];
    %Xk0=hformos(1,Zk,[1 0 1]);
    % Same thing
    Xk=abs(Hk(:)).^2./Sb(:);

    % Calculate residuals, removing values at the zero wavenumber and
    % wavenumbers above the minimum wavenumber "params.kiso"
    Xk(k==0)=NaN;
    Xk(k>params.kiso)=NaN;
    % An additional adjustment to try is reducing the annular radius of the
    % residuals a bit in case the data filter was not applied as quoted (recall 
    % conversation w Sofia on June 11 2025)
    %Xk(k>params.kiso.*0.8)=NaN;

    % And this should be the same thing again, except how it treats k=0
    [Lbar,~,~,momx,~,Xk1]=logliosl(k,thhat.*scl,params,Hk,1);
    Xk1=-Xk1-log(Sb(~~k));
    % The oldest way
    % Xkk1=-lkosl(k,thhat.*scl,params,Hk)-log(Sb);
    % difer(Xkk1(~~k)-Xk1,9,[],NaN)

    % Labeling  thing
    varibal='X';
    xstr2v=sprintf('quadratic residual 2%s',varibal);

    % Evaluate the likelihood
    % disp(sprintf('The loglihood is %8.3f',Lbar))

    %% TAKE CARE OF A FEW THINGS UPFRONT
    % Bounds of X residuals to show (functions as color axis in panel 3)
    boundsX0=[0 3*df];
    % Bounds of 2X residuals to how (functions as axis bounds in panel 1)
    bounds2X=2*boundsX0;
    % Color axis for the power spectral densities (panels 2 and 4)
    caxPSD=[-5 0]; % Venus only
    caxPSD=[-3.5 0]
    % Contours to be plotted on the power spectral density plots
    conPSD=[-5:-1];
    conPSD=[-3.5:-1] % Venus only
    % Get order of magnitude of last wavenumber for scaling
    om=round(log10(max(k(:))));

    % Set up the wavenumber/wavelength axes labels for panels 2-4
    xstr1=sprintf('x wavenumber (rad/%s) %s10^{%i}',unts,'\times',om);
    ystr1=sprintf('y wavenumber (rad/%s) %s10^{%i}',unts,'\times',om);
    xstr2=sprintf('x wavelength (%s)',unts);
    ystr2=sprintf('y wavelength (%s)',unts);

    % Prepare the overall figure
    fh=gcf;
    fig2print(fh,'landscape')
    set(fh,'PaperPositionMode','auto')

    % Select non-ridiculous values
    if any(isinf(Xk)); error('Adapt as in MLECHIPLOS and reconcile'); end
    allg=~isnan(Xk);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PANEL 1: TWICE CHI-SQUARED RESIDUALS AND QQ-PLOT [Top-left]
    axes(ah(1))

    % Get binning info for the twice chi-squared residuals histogram
    binos=5*round(log(length(Xk(allg))));
    % Rather this, to hit the grid lines
    binWidth=df/4;
    binos=(binWidth/2)*[1:2:bounds2X(2)*2/binWidth+1];
    [bdens,c]=hist(2*Xk(allg),binos);
    % Plot the histogram as a bar graph
    bdens=bdens/indeks(diff(c),1)/length(Xk(allg));
    bb=bar(c,bdens,1);
    % Each of the below should be df/2
    disp(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
		 varibal,nanmean(Xk(allg)),...
		 varibal,nanvar(Xk(allg))));

    set(bb,'FaceC',grey)
    hold on
    % Plot the ideal chi-squared distribution
    refs=linspace(0,max(bounds2X),100);
    plot(refs,chi2pdf(refs,df),'Linew',1.5,'Color','k')
    hold off

    % Labeling and cosmetic adjustments
    xlim(bounds2X)
    xl1(1)=xlabel(xstr2v);
    ylim([0 max(bdens)*1.05])
    yl1(1)=ylabel('probability density');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Prepare an overlay axis for the quantile-quantile plot
    ah2(1)=laxis(ah(1),0,0);
    axes(ah2(1))

    % Obtain qq-plot data
    % Note that the try/catch provides the necessary version upgrade
    try
        h=qqplot(2*Xk(allg),makedist('gamma','a',df/2,'b',2));
    catch
        h=qqplot(2*Xk,ProbDistUnivParam('gamma',[df/2 2]));
    end
    hx=get(h,'Xdata'); hx=hx{1};
    hy=get(h,'ydata'); hy=hy{1};
    delete(h)

    % Make the qq-plot
    qq0=plot(bounds2X,bounds2X,'k'); hold on
    qq=plot(hx,hy,'LineS','none','Marker','o','MarkerF','r',...
	    'MarkerE','r','MarkerS',2);

    % More cosmetic adjustments
    set(ah(1),'box','off')
    set(ah2(1),'xlim',get(ah(1),'xlim'),'yaxisl','r','box','off',...
       'xaxisl','t','color','none','ylim',get(ah(1),'xlim'))

    % Add labels and tickmarks
    ylr(1)=ylabel('quantile-quantile prediction');
    tkmks=bounds2X(1):2*df:bounds2X(2);
    set([ah(1) ah2(1)],'xtick',tkmks)
    set(ah2(1),'ytick',tkmks)

    % Add gridlines  % (c) Jos van der Geest
    try
        gh=gridxy([4 8],[4 8],'LineStyle',':');
    end

    % Information is power
    tbstr{1}=sprintf('%5.2f%%',100*sum(binWidth*bdens(1:end-1)));
    tbstr{2}=sprintf('max(2X) = %5.2f',max(2*Xk));
    tbstr{3}=sprintf('m(X) = %5.2f',nanmean(Xk));
    tbstr{4}=sprintf('v(X) = %5.2f',nanvar(Xk));

    % Calculate the mean of (X-df/2)^2 so it can be passed as output
    magx=nanmean((Xk(allg)-df/2).^2);
    %tbstr{5}=sprintf('mean([X-%d]^2) = %5.2f',df/2,magx);
    tbstr{5}=sprintf('s_X^2 = %5.2f',magx);

    % Do the test whether you accept this as a good fit, or not
    vr=8/sum(allg);
    [a,b,c,d]=normtest(magx,1,vr,0.05);
    % disp(sprintf('NORMTEST %i %5.3f %i',a,b,round(c)))
    if a==0; stp='accept'; else stp='reject'; end
    tbstr{6}=sprintf(' %s at %4.2f',stp,d);
    tbstr{7}=sprintf('p = %5.2f',b);

    % Give the x- and y-positions of the textbox
    tbx=repmat(bounds2X(2)-bounds2X(2)/30,1,length(tbstr));
    tby=[7.5 6.35-[0 1 2 3.5 4.5 5.5]];

    % Make the textbox(es) with unbordered fillboxes around them
    for i=1:length(tbstr)
        tb(i)=text(tbx(i),tby(i),tbstr{i},'HorizontalA','Right');
        gg=ext2lrtb(tb(i),1.1,1.1); delete(tb(i)); hold on
        fb(i)=fillbox(gg+[-0.1 0.08 0 0]);
        tb(i)=text(tbx(i),tby(i),tbstr{i},'HorizontalA','Right');
        hold off
        set(fb(i),'EdgeC','w','FaceC','w')
    end

    %%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PANEL 2: PREDICTED 2D POWER SPECTRUM [Top-right]
    axes(ah(2))

    % The top-left and bottom-right wavenumbers, scaled
    % The negatives are to trick IMAGEFNAN - we want a Nyquist, if even, at the top and the left,
    % and we want the zero right on the grid line
    % FJS Re-check this - the below does Venus right, but the old one (no negative, fliplr below)
    % looked right in the PDF of the fine samplings, but not on the screen - so
    % perhaps there is a rendering thing going on. Anyway, we are going with the
    % one where the small grids look right in the rendering, that is all
    c11=[kx(1)   -ky(1)]/10^om;
    cmn=[kx(end) -ky(end)]/10^om;

    % Remove the zero wavenumber value (to avoid issues with the
    % colorscale) and those at wavenumbers above the isotropic cutoff
    Sb(k==0)=NaN;
    Sb(k>params.kiso)=NaN;

    Sb=reshape(Sb,params.NyNx);

    % Scale to the largest predicted value and plot the power spectrum
    sclSb=log10(Sb./max(Sb(:)));
    imagefnan(c11,cmn,sclSb,'jet',caxPSD,[],[],0);

    % Label the wavenumber axes
    xl1(2)=xlabel(xstr1);
    yl1(2)=ylabel(ystr1);

    % For the wavelength axis, get tickmark values on the wavenumber axis
    xtk=get(ah(2),'xtick')*10^om;
    ytk=get(ah(2),'ytick')*10^om;

    if strcmp(unts,'km')
        rond=0;
        convo=1;
    elseif strcmp(unts,'mm')
        rond=2;
        convo=1;
    end

    % If in km, convert to m by dividing by 1000
    % Convert to wavelengths, recognizing the zero wavenumber
    xtkl=2*pi./xtk*convo;
    xtkl(isinf(xtkl))=[params.NyNx(2)-1]*params.dydx(2)*convo;
    ytkl=2*pi./ytk*convo;
    ytkl(isinf(ytkl))=[params.NyNx(1)-1]*params.dydx(1)*convo;

    % Create and label the wavelength axis
    [ah2(2),xl2(2),yl2(2)]=xtraxis(ah(2),xtk/10^om,round(xtkl,rond),xstr2,...
			           ytk/10^om,round(ytkl,rond),ystr2);

    % Return to the main axis and prepare to plot contours
    axes(ah(2)); hold on

    exx=xlim; eyy=ylim;
    
    % Calculate and plot contours
    %    [~,ch(1)]=contour(kx/10^om,fliplr(ky/10^om),sclSb,conPSD,'LineW',2);
    [cl{1},ch(1)]=contour(kx/10^om,ky/10^om,sclSb,conPSD,'LineW',2);
    cll{1}=clabel(cl{1});

    colormap(jet)
    caxis(caxPSD)

    axes(ah(2))
    axis([exx eyy])

    % Option to set the contours to black for visibility
    %set(hb,'color','k')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% PANEL 3: 2D X RESIDUALS [Bottom-left]
    axes(ah(3))

    % Plot the X residuals
    imagefnan(c11,cmn,reshape(Xk,params.NyNx),'gray',boundsX0,[],1,0);

    % The BOXTEX messes up the neat flat shading from imagefnan... or
    % does it?
    if ~isnan(params.kiso)
        % Place a box giving the wavelength of the isotropic wavenumber cutoff
        [~,zy]=boxtex('ll',ah(3),sprintf('%2.0f %s',2*pi./params.kiso,unts),...
		      12,[],[],1.1);
        set(zy,'FontSize',get(zy,'FontSize')-1)
    end

    % Label the wavenumber axes
    xl1(3)=xlabel(xstr1);
    yl1(3)=ylabel(ystr1);

    % Create and label the wavelength axis
    ah2(3)=xtraxis(ah(3),xtk/10^om,round(xtkl,rond),xstr2,...
	               ytk/10^om,round(ytkl,rond),ystr2);
    
    % Remove the right y-label to avoid clutter
    delete(get(ah2(3),'ylabel'))

    axes(ah(3))

    % Add a colorbar
    [cb,xcb]=addcb('vert',boundsX0,boundsX0,'gray',df,1);
    set(xcb,'string','quadratic residual X')
    moveh(cb,0.015)

    % Reposition the main axis
    set(ah(3),'position',[getpos(ah(3),1) getpos(ah(3),2) ...
		          getpos(ah(2),3) getpos(ah(2),4)])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PANEL 4: 2D TRUE POWER SPECTRUM [Bottom-right]
    axes(ah(4))

    % Calculate the periodogram of the edata
    Sk=abs(reshape(Hk(:),params.NyNx)).^2;

    % Remove the zero wavenumber value (to avoid issues with the
    % colorscale) and those at wavenubmers above the isotropic cutoff
    Sk(k==0)=NaN;
    Sk(k>params.kiso)=NaN;

    % Scale to the largest predicted value and plot the power spectrum
    sclSk=log10(Sk./max(Sb(:)));
    imagefnan(c11,cmn,sclSk,'jet',caxPSD,[],[],0);

    % Label the wavenumber axes
    xl1(4)=xlabel(xstr1);
    yl1(4)=ylabel(ystr1);

    % Create and label the wavelength axis
    ah2(4)=xtraxis(ah(4),xtk/10^om,round(xtkl,rond),xstr2,...
	           ytk/10^om,round(ytkl,rond),ystr2);

    % Return to the main axis and prepare to plot contours
    axes(ah(4)); hold on

    exx=xlim; eyy=ylim;
    
    % Place contours from the predicted power spectrum onto this true one
    [cl{2},ch(2)]=contour(kx/10^om,ky/10^om,sclSb,conPSD,'LineW',2);
    cll{2}=clabel(cl{2});
    caxis(caxPSD)

    axes(ah(4))
    axis([exx eyy])

    % Option to set the contours to black for visibility
    % set(hb,'color','k')

    %% FINAL COSMETIC ADJUSTMENTS

    % Adjust tickmarks
    longticks([ah ah2 cb])

    % Adjust the colorbar in panel 3
    axes(cb)
    set(cb,'YAxisLocation','right')
    moveh(cb,0.06)
    %  moveh(cb,0.01)
    cb.Position(2)=ah(3).Position(2);
    cb.Position(4)=ah(3).Position(4);
    cb.Position(3)=cb.Position(3)/1.75;
    % shrink(cb,1.3,1.27)

    % Adjust the main axes
    set(ah(2:4),'Box','On')

    % Give the overall figure a title
    axes(ah2(1))
    [~,b]=star69;
    if strcmp(b,'eggers8')
        spt=text(df,bounds2X(2)-df+1/2,stit);  
    else
        spt=title(ah2(1),stit);
        movev(spt,1)
    end

    movev([ah ah2 cb],-.02)
    movev([ah(1) ah(2) ah2(1) ah2(2)],0.01)

    % Set figure background color to white
    set(gcf,'color','W','InvertH','off')

    % Collect output
    vars={a,magx,ah,ah2,cb,ch,spt,cll};
    varargout=vars(1:nargout);
elseif strcmp(Hk,'demo1')
    % Image source: https://www.alexstrekeisen.it/english/meta/quartzite.php
    % Thin section, keep them inside a directory ROCKS and specify IFILES
    imnum=6; %2-7 are available
    warning off MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero
    fnam=fullfile(getenv('IFILES'),'ROCKS',sprintf('quartzite2012_%i.jpg',imnum));
    warning off MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero
    % Read the  image
    imdata=imread(fnam);
    immeta=imfinfo(fnam);
    % Initialize grid information
    p.NyNx=[immeta.Height immeta.Width];
    % can we use something more physical than pixels?
    unts='mm';
    % field of view [mm] of the image (height) where the width would have been 7mm
    fov=4.6638; 
    px=4.6638/immeta.Height;
    p.dydx=[px px];
    p.blurs=-1;

    % May need to decimate the image if it is too big
    dec=3;
    if any(p.NyNx>400)
        imdata=imdata(1:dec:end,1:dec:end,:);
        p.NyNx=[size(imdata,1) size(imdata,2)];
        p.dydx=p.dydx*dec;
    end
    % Convert to gray scale
    try
        Hx=double(im2gray(imdata));
    catch
        Hx=[double(imdata(:,:,1))+double(imdata(:,:,2))+double(imdata(:,:,3))]/3;
    end
    % Remove constant
    km=mean(Hx(:));
    Hx=Hx(:)-km;

    % Set up the call for the main MLECHIPSDOSL routine, including applying the
    % smooth taper for the estimation
    p.blurs=-1; p.kiso=NaN; 
    pt=p; Tx=gettaper(pt,'cosine',0.10);
    pt.taper=Tx;
    % Make the estimate and request covariance using the sampling method
    % Don't redo if you had it
    fname=fullfile(getenv('IFILES'),'HASHES',sprintf('%s-%s',upper(mfilename),hash(Hx,'SHA-1')));
    if ~exist(sprintf('%s.mat',fname),'file')
        % You could give it something close either in thini or aguess if you've run it before
        % but this doesn't always actually make it easier
        thini=[];
        thhat=NaN;
        while isnan(thhat)
            [thhat,covFHhJ,~,scl,~,~,Hk]=mleosl(Hx,thini,pt,[],[],[],[],[],1);
        end
        save(fname,'thhat','scl','Hk','covFHhJ')
    else
        load(fname)
    end
        
    % Best fits for a certain image and decimation
    % im dec
    % 2 2 4043.5        1.2396      0.031736
    % 3 2 3848.7        1.3669      0.031397
    % 4 2 4172.9        1.2024      0.035007
    % 5 2 4425.0        1.2230      0.035612
    % 6 2
    % 7 2

    % 2 3 4823.7        0.93916     0.045664 
    % 3 3 4458.1        1.07860     0.041598
    % 4 3 4928.9        0.92742     0.049980
    % 5 3 5250.8        0.93905     0.050766
    % 6 3 4876.3        0.90899     0.051175  ACCEPT
    % std  196          0.0094      0.001532
    % 7 3 5237.7        0.93482     0.050461
    
    % Simulate at the estimate, without the taper
    p.blurs=Inf;
    HxS=simulosl(thhat.*scl,p);

    % Make a first figure with the oboserved and the simulated field
    f1=figure(1);
    fig2print(f1,'landscape')
    set(f1,'PaperPositionMode','auto')
    tts={'data','synthetic'};
    unts='mm';
    axs={sprintf('dimension 1 (%s)',unts) 
         sprintf('dimension 2 (%s)',unts)};
    fw='normal';
    % Common color range
    cmap=gray;
    cax=prctile([Hx ; HxS],[5 95])

    clf
    [ah,ha,H]=krijetem(subnum(2,2));
    % Plot the data
    [tl(1),xl(1),yl(1),cb1(1)]=plotit(flipud(v2s(Hx,p)),p,ah(1),tts{1},axs,cmap,cax,fw,km);
    % Plot the synthetic
    [tl(2),xl(2),yl(2),cb1(2)]=plotit(       v2s(HxS,p),p,ah(3),tts{2},axs,cmap,cax,fw,km);

    label(ah([1 3 2 4]),'ul',[],0,0)
    
    % Just prepare for LaTeX trim and clip on 8.5.0.197613 (R2015a)
    delete(ah([2 4]))

    % Labels etc
    cb1(1).XLabel.String='8-bit grey scale';
    cb1(2).XLabel.String='dd';

    delete(cb1(2))
    delete(tl([1 2]))
    allofem=findobj('FontSize',10);
    for ondex=1:length(allofem)
        set(allofem(ondex),'FontSize',9)
    end
    
    figure(2)
    clf
    % Remember that MLEOSL already returned variance scaled data
    scl2=scl; scl2(1)=1;
    [~,~,nah,nah1,cb,ch,spt,cll]=mlechipsdosl(Hk,thhat,scl2,pt,...
                         sprintf('%s = %i %s = %4.2f  %s = %4.2f %s',...
                        '\sigma^2',round(thhat(1)*scl(1)),...
                        '\nu',thhat(2)*scl(2),...
                        '\rho',round(thhat(3)*scl(3),3),unts),...
                         [],unts);
    [bh,th]=label([nah1(1) nah([3 2 4])],'ul',[],2,0);

    % Now prettify the contour labels
    % blablabl
    % we'll grab 'xdata' and 'ydata'
    % fillbox(ext2lrtb(ext2lrtb(cll{2}(4)))
keyboard
    
    figure(1)
    figna=figdisp([],sprintf('%s_%i','demo1_1',imnum),[],1);
    figure(2)
    figna=figdisp([],sprintf('%s_%i','demo1_2',imnum),[],1);
elseif strcmp(Hk,'demo2')
    % Bathymetry from GEBCO
    % Data window selected as a subset of the region in DEMO2
    fnam=fullfile(getenv('IFILES'),'GEBCO','GEBCO2024Grid_AtlanticMERMAIDpatch.mat');
    dat=load(fnam);
    % compare
    [LOLO,LALA]=meshgrid(dat.mlons,dat.mlats);
    % Note that this will triggger parallel pool
    % z2008=gebco(LOLO,LALA,2008);
    % z2014=gebco(LOLO,LALA,2014);
    z2019=gebco(LOLO,LALA,2019);
    % plot(z2019(:),dat.sgelev(:),'.')
    % Decide to keep old GEBCO
    % elev=dat.sgelev;
    elev=z2019;
    lat=dat.mlats;
    lon=dat.mlons;
    p=dat.p;
    % Crop the region?
    ndx=500:numel(lon)-1.2e3;
    tdx=300:1000;
    imnum=1;
    % May need to decimate the field if it is too big
    if any(p.NyNx>400)
        dec=3;
        tdx=tdx(1:dec:end);
        ndx=ndx(1:dec:end);
        p.dydx=p.dydx.*dec;
    end
    p.NyNx=[numel(tdx) numel(ndx)];
    p.lon=lon;
    p.lat=lat;
    p.ndx=ndx;
    p.tdx=tdx;
    % Possibly no need for further conversion
    Hx=double(elev(tdx,ndx));
    
    % Remove regional trend
    [~,~,~,~,~,Z]=planefit(v2s(Hx,p));
    km=mean(Hx(:));
    Hx=Hx(:)-Z(:);

    % Set up the call for the main MLECHIPSDOSL routine, including applying the
    % smooth taper
    p.blurs=-1;
    p.kiso=NaN;
    pt=p;
    Tx=gettaper(pt,'cosine',0.10);
    pt.taper=Tx;
    % Make the estimate and request covariance using the sampling method
    % Don't redo if you had it 
    fname=fullfile(getenv('IFILES'),'HASHES',sprintf('%s-%s',upper(mfilename),hash(Hx,'SHA-1')));
    if ~exist(sprintf('%s.mat',fname),'file')
        % You could give it something close either in thini or aguess if you've run it before
        % but this doesn't always actually make it easier
        thini=[];
        thhat=NaN;
        while isnan(thhat)
            [thhat,covFHhJ,~,scl,~,~,Hk]=mleosl(Hx,[],pt,[],[],[],[],[],1);
        end
        save(fname,'thhat','scl','Hk','covFHhJ')
    else
        load(fname)
    end

    % With deplaning dec 3
    % Estimated theta :  46139        1.2558        2.4795
    % FishJFish. std :      1.25e+03       0.00782      0.038236

    % Simulate at the estimate, without the taper
    p.blurs=Inf;
    HxS=simulosl(thhat.*scl,p);

    % Make a first figure with the observed and the simulated field
    f1=figure(1);
    fig2print(f1,'landscape')
    set(f1,'PaperPositionMode','auto')
    tts={'data','synthetic'};
    axs={sprintf('longitude (%s)',str2mat(176)) 
         sprintf('latitude (%s)',str2mat(176))};
    fw='normal';
    % Common color range
    [cmap,~,ziro]=sergeicol;
    cmap=cmap(1:ziro,:);
    cax=prctile([Hx ; HxS],[0 100])+km;
    % Lighten up by subtracting right
    cax=cax+[0 -500];

    clf
    [ah,ha,H]=krijetem(subnum(2,2));
    % Plot the data but keep the mean again
    [tl(1),xl(1),yl(1),cb1(1)]=plotit2(flipud(v2s(Hx,p))+km,p,ah(1),tts{1},axs,cmap,cax,fw);
    % Plot the synthetic
    [tl(2),xl(2),yl(2),cb1(2)]=plotit2(       v2s(HxS,p)+km,p,ah(3),tts{2},axs,cmap,cax,fw);

    label(ah([1 3 2 4]),'ul',[],0,0)

    % Just prepare for LaTeX trim and clip on 8.5.0.197613 (R2015a)
    delete(ah([2 4]))

    % Labels etc
    cb1(1).XLabel.String='8-bit grey scale';
    cb1(2).XLabel.String='dd';

    delete(cb1(2))
    delete(tl([1 2]))
    allofem=findobj('FontSize',10);
    for ondex=1:length(allofem)
        set(allofem(ondex),'FontSize',9)
    end
    
    figure(2)
    clf
    scl2=scl;scl2(1)=1;
    [~,~,nah,nah1,cb,ch,spt,cll]=mlechipsdosl(Hk,thhat,scl2,pt,...
                         sprintf('%s = %i %s = %4.2f  %s = %4.2f %s',...
                        '\sigma^2',round(thhat(1)*scl(1)),...
                        '\nu',thhat(2)*scl(2),...
                        '\rho',round(thhat(3)*scl(3),3),unts),...
                         [],unts);

    [bh,th]=label([nah1(1) nah([3 2 4])],'ul',[],2,0);

    % Now prettify the contour labels
    % blablabl
    % we'll grab 'xdata' and 'ydata'
    % fillbox(ext2lrtb(ext2lrtb(cll{2}(4)))

    keyboard

    figure(1)
    figna=figdisp([],sprintf('%s_%i','demo2_1',imnum),[],1);
    figure(2)
    figna=figdisp([],sprintf('%s_%i','demo2_2',imnum),[],1);
elseif strcmp(Hk,'demo3')
    % Fault roughness examples from data at https://doi.org/10.5281/zenodo.6411819
    % associated with 10.1029/2022JB025511
    % parsed and reformatted using readfault.m as mat files
    % ddir='~/Downloads/guerin_et_al_2022_GRL_data/roughness_analysis'; 
    fnams={'R1_before_bot','R1_after_bot',...
           'R1_before_top','R1_after_top',...
           'S2_before_bot','S2_after_bot',...
           'S2_before_top','S2_after_top'};
    foril=1;
    if ~foril
        % Load for exploratory purposes
        for ind=1:numel(fnams);
            tmp=load(fullfile(getenv('IFILES'),'TRIBOLOGY',sprintf('%s%s%s','olw_',fnams{ind},'.mat')));
            block{ind}=tmp.dat;
            p{ind}=[];
            p{ind}.NyNx=size(block{ind});
        end

        % Quick visual comparison of the eight block samples
        clf; [ah,ha,H]=krijetem(subnum(numel(fnams)/2,2));
        for ind=1:numel(fnams);
            axes(ah(ind));
            imagefnan([1 p{ind}.NyNx(2)],[1 p{ind}.NyNx(1)],...
                      block{ind}-nanmean(block{ind}),[],[-0.01 0.01])
            title(strrep(sprintf('%s',fnams{ind}),'_',' '))
        end

        % The smooth experiments taken as full rectangular observations or as
        % patches estimate values of rho that are smaller than the sample spacing.
        % Spatially, the smooth surfaces possess many blunt peaks and deep pits,
        % while the roughened surfaces strike a better balance, though there is more
        % structure that could possibly be modeled by a higher order polynomial.
        % All samples possess high gradients between extreme peaks and troughs that
        % are a challenge to capture in the modeled field.
    else
        % Load for real
        % Choose which block to work with:
        fnum=3; 
        % Hx=block{fnum};
        tmp=load(fullfile(getenv('IFILES'),'TRIBOLOGY',sprintf('%s%s%s','olw_',fnams{fnum},'.mat')));
        Hx=tmp.dat;
    end

    % Set the grid's parameters so that the grid size is 100mm by 50mm
    p=[];
    p.NyNx=size(Hx);
    p.dydx=[50 100]./p.NyNx;
    unts='mm';

    % The data itself is measured from an elliptical cross-section, so we will
    % first crop the dataset to have a rectangular field of observation
    decm=1;
    cy=225:decm:p.NyNx(1)-225;
    cx=450:decm:p.NyNx(2)-450;

    % Now we might choose to work with the full field, use default values for a
    % few of the datasets, or explore the dataset and capture a patch of
    % interest based on the displayed axes limits (explr==1)

    % Option to pause and look around
    explr=0;
    if fnum==1 & ~explr
      xl=[1163 1348];yl=[80   250];
    elseif fnum==3 & ~explr
      xl=[-30 76]; yl=[406 513];
    else
      disp('zoom in to select a data patch for analysis')
      figure(3); clf
      hi=imagesc(1:numel(cy),1:numel(cx),Hx(cy,cx));
      colorbar
      h = zoom();
      hpan = pan(gcf);
      h.ActionPostCallback = @changecolorbar;
      hpan.ActionPostCallback = @changecolorbar;
      keyboard
      yl=floor(xlim);xl=floor(ylim);
    end

    % Adjust the crop
    cy=cy(1)+yl(1):cy(1)+yl(2);
    cx=cx(1)+xl(1):cx(1)+xl(2);
    % Take the cropped data and adjust the grid parameters 
    Hx=Hx(cy,cx);
    p.NyNx=size(Hx)
    p.dydx=p.dydx*decm;

    % Remove a first order polynomial
    [~,~,~,~,~,Z]=planefit(v2s(Hx,p));
    Hx=Hx(:)-Z(:);
    Hx=Hx(:)-nanmean(Hx(:));

    % Store the dataset as Hx0
    Hx0=Hx;
    % Fill the extreme peaks and troughs outside of mnprc -- mxprc in the
    % dataset that will be used for the MLE analysis
    mnprc=7.5;
    mxprc=92.5;
    Hx(Hx<prctile(Hx,mnprc))=prctile(Hx,mnprc);
    Hx(Hx>prctile(Hx,mxprc))=prctile(Hx,mxprc);

    % Set up the call for the main MLECHIPSDOSL routine, including applying the
    % smooth taper
    p.blurs=-1; p.kiso=NaN; 
    % Setting kiso makes the residual pattern less complicated, but it does not
    % change the estimate: pi/max(p.dydx);
    pt=p; Tx=gettaper(pt,'cosine',0.10);
    pt.taper=Tx;
    % Make the estimate and request covariance using the sampling method
    % Don't redo if you had it
    fname=fullfile(getenv('IFILES'),'HASHES',sprintf('%s-%s',upper(mfilename),hash(Hx,'SHA-1')));
    if ~exist(sprintf('%s.mat',fname),'file')
        % You could give it something close either in thini or aguess if you've run it before
        % but this doesn't always actually make it easier
        thini=[];
        thhat=NaN;
        while isnan(thhat)
            [thhat,covFHhJ,~,scl,~,~,Hk]=mleosl(Hx,thini,pt,[],[],[],[],[],1);
        end
        save(fname,'thhat','scl','Hk','covFHhJ')
    else
        load(fname)
    end
    
    % Simulate at the estimate, without the taper
    p.blurs=Inf;
    HxS=simulosl(thhat.*scl,p);

    % We will be setting the color axes for the data and the synthetic. 
    % A few options to pick from:
    caxmth=3;
    switch caxmth
      case 1
        % halverange based; this is symmetric, which we may not want here
        % share the same caxis determined from the truncated data
        cax =halverange(Hx,95);
        caxS=cax;
      case 2
        % percentile based
        % share the same caxis determined from the truncated data
        cax =prctile(Hx,[5 95]);
        caxS=cax;
     case 3 
        % percentile based
        % set distinctly for the truncated data and the synthetic
        cax =prctile(Hx ,[5 95]);
        caxS=prctile(HxS,[5 95]);
    end

    % Make a first figure with the oboserved and the simulated field
    f1=figure(1);
    fig2print(f1,'landscape')
    set(f1,'PaperPositionMode','auto')
    tts={'data','synthetic'};
    unts='mm';
    axs={sprintf('dimension 1 (%s)',unts)
         sprintf('dimension 2 (%s)',unts)};
    fw='normal';
    cmap=pink;
    cmap=cmap(2:end-10,:);

    clf
    [ah,ha,H]=krijetem(subnum(2,2));
    % Plot the data
    [tl(1),xl(1),yl(1),cb1(1)]=plotit(flipud(v2s(Hx,p)),p,ah(1),tts{1},axs,cmap,cax,fw);
    % Plot the synthetic
    [tl(2),xl(2),yl(2),cb1(2)]=plotit(       v2s(HxS,p),p,ah(3),tts{2},axs,cmap,caxS,fw);

    label(ah([1 3 2 4]),'ul',[],0,0)

    % Just prepare for LaTeX trim and clip on 8.5.0.197613 (R2015a)
    delete(ah([2 4]))

    % Run fnum 3 for the paper figures
    % True theta :    4.3413e-05        1.6389       0.10333
    %                 4.13e-06          0.0228       0.0042244

    % Labels etc
    cb1(1).XLabel.String='8-bit grey scale';
    cb1(2).XLabel.String='dd';

    delete(cb1(2))
    delete(tl([1 2]))
    allofem=findobj('FontSize',10);
    for ondex=1:length(allofem)
        set(allofem(ondex),'FontSize',9)
    end
    
    figure(2)
    clf
    % Remember that MLEOSL already returned variance scaled data
    scl2=scl; scl2(1)=1;
    [~,~,nah,nah1,cb,ch,spt,cll]=mlechipsdosl(Hk,thhat,scl2,pt,...
                         sprintf('%s = %i %s = %4.2f  %s = %4.2f %s',...
                        '\sigma^2',round(thhat(1)*scl(1)),...
                        '\nu',thhat(2)*scl(2),...
                        '\rho',round(thhat(3)*scl(3),3),unts),...
                         [],unts);
    [bh,th]=label([nah1(1) nah([3 2 4])],'ul',[],2,0);

    % Now prettify the contour labels
    % blablabl
    % we'll grab 'xdata' and 'ydata'
    % fillbox(ext2lrtb(ext2lrtb(cll{2}(4)))

    keyboard
    figure(1)
    figna=figdisp([],sprintf('%s_%i','demo3_1',fnum),[],1);
    figure(2)
    figna=figdisp([],sprintf('%s_%i','demo3_2',fnum),[],1);

    % Before the experiment, for the rough block cropped to its center, no deplane:
    % fnum 1, th [0.000118        0.754      106.38],  FishJFish. std : [2.63e-06     0.01        2.446]
    % fnum 3, th [0.0001381       0.6474     117.78],  FishJFish. std : [3.4e-06      0.00845     3.2111]
    % Before the experiment, smooth: 
    % fnum 5, th [ 9.0964e-06     0.68947     36.153], FishJFish. std : [ 7.37e-08    0.0275      0.50783]
    % fnum 7, th [ 9.6435e-06     0.72163     40.372], FishJFish. std : [ 9.02e-08    0.0245      0.6033]
    % After the experiment, rough:
    % fnum 2, th [0.00013773      0.61311    118.61],  FishJFish. std : [ 3.22e-06    0.00837     3.2071]
    % fnum 4, th [0.00011017      0.55212    154.12],  FishJFish. std : [ 3.21e-06    0.00714     5.4017]
    % After the experiment, smooth:
    % fnum 6, th [  9.75e-06      0.557       42.951], FishJFish. std : [ 9.51e-08    0.0177      0.74717]
    % fnum 8, th [ 1.1422e-05     0.6619      43.266], FishJFish. std : [ 1.13e-07    0.0211      0.68681]
elseif strcmp(Hk,'demo4')
    % Load a patch of Venus; the ``secret'' to passing the test at this point is
    % decimating the data, or using blurs>-1. We don't want to do the latter, so
    % choose a patch that is large to begin with. Below are some options:
    imnum=2;
    % pth='~/Documents/Venus/';
    pth=fullfile(getenv('IFILES'),'VENUS','DATA','plmData');
    fnam='plmVenus_D-5.mat';
    if imnum==1
      topo=modload(fullfile(pth,fnam),'V0010_03');
      tdx1=1;      ndx1=1;
      [tdx2,ndx2]=size(topo.dataP.dp);
    elseif imnum==2
      topo=modload(fullfile(pth,fnam),'V0022_03');
      tdx1=1;      ndx1=1;
      [tdx2,ndx2]=size(topo.dataP.dp);
    elseif imnum==3
      topo=modload(fullfile(pth,fnam),'V0066_03');
      tdx1=1;      ndx1=1;
      [tdx2,ndx2]=size(topo.dataP.dp);
    end
    unts=topo.dataP.unit;
    p=topo.params;

    % It would be nice to work with km and not m; make the conversion in dydx
    unts='km';
    p.dydx=p.dydx/1e3;
    elev=topo.dataP.dp;
    % Decimate
    dec=2;
    p.dydx=p.dydx.*dec;
    tdx=tdx1:dec:tdx2;
    ndx=ndx1:dec:ndx2;
    Hx=elev(tdx,ndx);
    p.NyNx=size(Hx);

    % Remove any planar trends
    [~,~,~,~,~,Z]=planefit(v2s(Hx,p));
    km=mean(Hx(:));
    Hx=Hx(:)-Z(:);
    Hx=Hx-mean(Hx);
    Hx=Hx(:);

    % Option to fill the extreme peaks and troughs outside of mnprc -- mxprc in
    % the dataset that will be used for the MLE analysis
    % This is not so important to these datasets
    mnprc=1;
    mxprc=99;
    Hx(Hx<prctile(Hx,mnprc))=prctile(Hx,mnprc);
    Hx(Hx>prctile(Hx,mxprc))=prctile(Hx,mxprc);
    Hx=Hx-mean(Hx);

    % Set up the call for the main MLECHIPSDOSL routine, including applying the
    % smooth taper
    p.blurs=-1; 
    p.kiso=pi./max(p.dydx);
    pt=p; Tx=gettaper(pt,'cosine',0.1);
    pt.taper=Tx;
    % Make the estimate and request covariance using the sampling method
    % Don't redo if you had it
    fname=fullfile(getenv('IFILES'),'HASHES',sprintf('%s-%s',upper(mfilename),hash(Hx,'SHA-1')));
    if ~exist(sprintf('%s.mat',fname),'file')
        % You could give it something close either in thini or aguess if you've run it before
        % but this doesn't always actually make it easier
        thini=[];
        thhat=NaN;
        while isnan(thhat)
            [thhat,covFHhJ,~,scl,~,~,Hk]=mleosl(Hx,[],pt,[],[],[],[],[],1);
        end
        save(fname,'thhat','scl','Hk','covFHhJ')
    else
        load(fname)
    end
    % Simulate at the estimate, without the taper
    p.blurs=Inf;
    HxS=simulosl(thhat.*scl,p);

    % Make a first figure with the oboserved and the simulated field
    f1=figure(1);
    fig2print(f1,'landscape')
    set(f1,'PaperPositionMode','auto')
    tts={'data','synthetic'};
    axs={sprintf('transformed longitude (%s)',str2mat(176))
         sprintf('transformed latitude (%s)',str2mat(176))};
    fw='normal';
    % Common color range
    cmap=flipud(kelicol);
    % Be like Gabe? With the planefit, his 5 95 is too restrictive
    perx=[2.5 97.5];
    perx=[2.5 97.5];
    perx=[1 99];
    cax1=prctile(Hx ,perx);
    cax1=prctile([Hx HxS],perx);
    cax2=prctile(HxS,perx); cax2=cax1;

    clf
    [ah,ha,H]=krijetem(subnum(2,2));

    p.lon=topo.geo.lonrDx; p.lat=topo.geo.latrDx; p.ndx=[1 2]; p.tdx=[1 2];
    % Plot the data but keep the mean again
    [tl(1),xl(1),yl(1),cb1(1)]=plotit2(flipud(v2s(Hx,p))+km,p,ah(1),tts{1},axs,cmap,cax1,fw);
    % Plot the synthetic
    [tl(2),xl(2),yl(2),cb1(2)]=plotit2(       v2s(HxS,p)+km,p,ah(3),tts{2},axs,cmap,cax2,fw);

    label(ah([1 3 2 4]),'ul',[],0,0)

    % Just prepare for LaTeX trim and clip on 8.5.0.197613 (R2015a)
    delete(ah([2 4]))

    % Labels etc
    cb1(1).XLabel.String='8-bit grey scale';
    cb1(2).XLabel.String='dd';

    delete(cb1(2))
    delete(tl([1 2]))
    allofem=findobj('FontSize',10);
    for ondex=1:length(allofem)
        set(allofem(ondex),'FontSize',9)
    end

    figure(2)
    clf
    scl2=scl;scl2(1)=1;
    [~,~,nah,nah1,cb,ch,spt,cll]=mlechipsdosl(Hk,thhat,scl2,pt,...
                         sprintf('%s = %i %s = %4.2f  %s = %4.2f %s',...
                        '\sigma^2',round(thhat(1)*scl(1)),...
                        '\nu',thhat(2)*scl(2),...
                        '\rho',round(thhat(3)*scl(3),3),unts),...
                         [],unts);
    [bh,th]=label([nah1(1) nah([3 2 4])],'ul',[],2,0);

    % Now prettify the contour labels
    % blablabl
    % we'll grab 'xdata' and 'ydata'
    % fillbox(ext2lrtb(ext2lrtb(cll{2}(4)))

    keyboard

    figure(1)
    figna=figdisp([],sprintf('%s_%i','demo4_1',imnum),[],1);
    figure(2)
    figna=figdisp([],sprintf('%s_%i','demo4_2',imnum),[],1);
elseif strcmp(Hk,'demo5')
    % WILL REQUIRE CLEANUP WITH CLABELS ETC
    % Sea surface height anomaly, near real time
    % Source: https://doi.org/10.48670/moi-00149
    % ``Altimeter satellite gridded Sea Level Anomalies (SLA) computed with respect to a
    % twenty-year [1993, 2012] mean. The SLA is estimated by Optimal Interpolation,
    % merging the L3 along-track measurement from the different altimeter missions
    % available. Part of the processing is fitted to the Global Ocean.''
    %
    % Level 4 processing -> filtering? see this document: 
    % https://documentation.marine.copernicus.eu/QUID/CMEMS-SL-QUID-008-032-068.pdf

    % Data window selected to correspond with proposed Atlantic MERMAID positions 
    pth='~/Documents/OceanFloor/SSHA/Data/';
    fnam=append(pth,...
                'cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1750366404080.nc');
    sla=ncread(fnam,'sla')';
    lat=ncread(fnam,'latitude');
    lon=ncread(fnam,'longitude');
    % Needs to be at least within these indices to avoid land
    ndx=40:numel(lon)-177;
    tdx=1:numel(lat);
    % But, we might want something a bit more square for making a nice visual,
    % so trim it back further 
    ndx=numel(lon)/2:numel(lon)-177;
    % And again to get away from spatially varying means
    tdx=1:numel(lat)*.63;
    % Better?
    Hx=sla(tdx,ndx);
    Hx=Hx(:)-mean(Hx(:));
    p=[];
    p.NyNx=[numel(tdx) numel(ndx)];
    % units in degrees
    p.dydx=[unique(lat(tdx(2):tdx(end))-lat(tdx(1):tdx(end-1)))...
            unique(lon(ndx(2):ndx(end))-lon(ndx(1):ndx(end-1)))];
    p.dydx=double(p.dydx);
    % units in km
    p.dydx=deg2km(p.dydx);
    p.blurs=-1;
    %
    figure(1);clf
    [ah,ha,H]=krijetem(subnum(2,3));
    axes(ah(1))
    [~,cax]=imagefnan([lon(ndx(1)) lat(tdx(1))],[lon(ndx(end)) lat(tdx(end))],...
                      v2s(Hx,p),'sky')
    xlabel(sprintf('longitude (%s)',char(176)))
    ylabel(sprintf('latitude (%s)',char(176)))
    ti(1)=title('Sea Surface Height Anomaly','FontWeight','normal');
    longticks
    % Set up the call for the main MLECHIPSDOSL routine, including applying the
    % smooth taper
    % TODO: set kiso?
    p.blurs=-1; p.kiso=0.5*pi./max(p.dydx);%NaN;
    pt=p;
    Tx=gettaper(pt,'cosine',0.10);
    pt.taper=Tx;
    thhat=NaN; while isnan(thhat); [thhat,~,~,scl,~,~,Hk,k]=mleosl(Hx,[],pt); end
    th=thhat.*scl;
    p.blurs=Inf;
    th=thhat.*scl;
    HxS=simulosl(th,p);
    tlabs=append('$\sigma=$ ',sprintf('%0.2f',sqrt(th(1))),...
                 ', $\nu=$ ',sprintf('%0.1f',th(2)),...
                 ', $\rho=$ ',sprintf('%0.2f km',th(3)),' $|$ ',...
                 sprintf('%0.1f',p.NyNx(1)*p.dydx(1)),'km x ',...
                 sprintf('%0.1f',p.NyNx(2).*p.dydx(2)),'km');
    % about [0.0014    7.9832   32.5332] in units [m^2 -- km] without taper
    % about [0.0022    3.5045   51.6191] in units [m^2 -- km] with smooth taper
    % looks good:
    p.blurs=Inf;
    axes(ah(4))
    imagefnan([lon(ndx(1)) lat(tdx(1))],[lon(ndx(end)) lat(tdx(end))],...
              v2s(HxS,p),'sky',cax)
    xlabel(sprintf('longitude (%s)',char(176)))
    ylabel(sprintf('latitude (%s)',char(176)))
    longticks
    ti(4)=title('Synthetic','FontWeight','normal');
    p.blurs=-1;
    figure(2); clf
    scl2=scl;scl2(1)=1;
    [~,~,nah,nah1,cb,ch,spt]=mlechipsdosl(Hk,thhat,scl2,pt,tlabs);
    for ind=1:4
        nnd=mod(ind,3)+floor(ind/3)*4+1;
        cah(ind)=copyobj(nah(ind),figure(1));
        set(cah(ind),'Position',get(ah(nnd),'position'))
        if ind==3
            axes(cah(ind))
            ncb=colorbar();
            ncb.Location='eastoutside';
            set(cah(ind),'Position',get(ah(nnd),'position'))
            ccb=copyobj(cb,figure(1));
            set(ccb,'Position',get(ncb,'position'))
            delete(ncb)
            moveh(ccb,0.02)
        elseif ind==2|ind==4
            colormap(jet)
        end
        cah1(ind)=copyobj(nah1(ind),figure(1));
        set(cah1(ind),'Position',get(ah(nnd),'position'))
        delete(ah(nnd))
    end
    f1ch=get(figure(1),'Children');
    set(f1ch(8).Title,'Interpreter','latex')
    % TODO: CORRECT THE WAVENUMBER AXES 
    for xnd=[2 3 4]
        cah1(xnd).XLabel.String='x wavelength (km)';
        xtkl=2*pi./get(cah1(xnd),'xtick');
        xtkl(isinf(xtkl))=[p.NyNx(2)-1]*p.dydx(2);
        ytkl=2*pi./get(cah1(xnd),'xtick');
        ytkl(isinf(ytkl))=[p.NyNx(1)-1]*p.dydx(1);
        cah1(xnd).XTickLabel=round(xtkl,1);
        cah1(xnd).YTick=cah1(xnd).XTick;
        cah1(xnd).YTickLabel=round(ytkl,1);
        if xnd~=3; cah1(xnd).YLabel.String='y wavelength (km)'; end
    end
    keyboard
    movev([ti(1) ti(4)],0.12)
    figna=figdisp([],'demo5_kiso',[],1);
end

% Just the space plots
function [tl,xl,yl,cb]=plotit(d,p,ah,tts,axs,cmap,cax,fw,km)
% Must make this active
axes(ah)
imagefnan([1 p.NyNx(1)],[p.NyNx(2) 1],d,cmap,cax);
axis image
longticks(ah)
ah.XTick=     [1 p.NyNx(2)]+[-1 1]*0.5;
ah.XTickLabel=round([0 p.NyNx(2)].*p.dydx(2));
ah.YTick=[1 p.NyNx(1)]+[-1 1]*0.5;
ah.YTickLabel=round([0 p.NyNx(1)].*p.dydx(1));
tl=title(tts,'FontWeight',fw);
xl=xlabel(axs{1});
yl=ylabel(axs{2});
movev(tl,range(ylim)/20);
% Put a colorbar
goodpos=[getpos(ah,1) getpos(ah,2)-getpos(ah,4)/5 getpos(ah,3) getpos(ah,4)/15];

% cb=colorbarf('hor',10,'Helvetica',goodpos);
% If it's funky FLOOR the first one, CEIL the last one, and ROUND the rest
cb=addcb(goodpos,cax,cax,cmap,round(km+linspace(min(cax),max(cax),6)));

% Just the space plots in untethered space coordinates
function [tl,xl,yl]=plotit2(d,p,ah,tts,axs,cmap,cax,fw)
% Must make this active
axes(ah)
imagefnan([p.lon(p.ndx(1)) p.lat(p.tdx(1))],[p.lon(p.ndx(end)) p.lat(p.tdx(end))],d,cmap,cax)
axis image
longticks(ah)
ext=[p.lon(p.ndx(1)) p.lon(p.ndx(end))];
ah.XTick=ext;
ah.XTickLabel=round(ext,2);
wait=[p.lat(p.tdx(1)) p.lat(p.tdx(end))];
ah.YTick=wait;
ah.YTickLabel=round(wait,2);
tl=title(tts,'FontWeight',fw);
xl=xlabel(axs{1});
yl=ylabel(axs{2});
movev(tl,range(ylim)/20);

% For exploring the demo4 datasets
function changecolorbar(src, ~)
% straight from https://www.mathworks.com/matlabcentral/answers/1449774-dynamic-colorbar-change-with-window-size-corresponding-with-different-data-area
h = gco;
xl = xlim; 
yl = ylim;
x = find(h.XData>=xl(1) &  h.XData<=xl(2));
iy = find(h.YData>=yl(1) &  h.YData<=yl(2));
C = h.CData(ix,iy);
caxis([min(C(:)) max(C(:))]);
