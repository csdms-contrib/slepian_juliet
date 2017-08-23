function varargout=mlechipsdosl(Hk,thhat,params,stit,ah)
% [mag,ah,ah2,cb,ch,spt]=MLESCHIPSDOSL(Hk,thhat,params,stit,ah)
%
% Makes a plot of the quadratic residuals and their interpretation for a
% Matern likelihood model of a single-field Fourier-domain data patch, as
% well as plots of the actual power spectral density of that data patch and
% its predicted power spectral density based on the model. These residuals
% should have a chi-squared distribution, so the produced figure is useful
% in determining how well the Matern likelihood model fits the data.
% 
% The figure produced consists of four panels:  
%    (TL) Histogram and qq-plot of the chi-squared residuals
%    (TR) Predicted 2D power spectral density based on "thhat"
%    (BL) 2D chi-squared residuals
%    (BR) "True" power spectral density with contours from (TR) overlain
%
% INPUT:
%
% Hk         The Fourier-domain data, e.g. from SIMULOSL
% thhat      The three Matern parameters
% params     Parameter set pertaining to the data, e.g. from SIMULOSL
% stit       Title for the overall figure [defaulted]
% ah         A 4x4 group of axis handles [defaulted]
%
% OUTPUT:
%
% mag         The "magic" parameter which I've come to call s_X^2
% ah,ah1      Main and secondary axis (panels 1-4) handles
% cb          Colorbar (panel 3) handle
% ch          Contour (panels 2 and 4) handles
% spt         Supertitle handle
%
% EXAMPLES:
%
% mlechipsdosl('demo1') % Runs itself for an example and a picture
%
% See also BLUROS, IMAGEFNAN, MATERNOS, LOGLIOSL, SIMULOSL, QQPLOT
%
% Last modified by gleggers-at-princeton.edu, 04/17/2014
% Last modified by fjsimons-at-alum.mit.edu, 07/08/2015

% Default values for conducting demos
defval('Hk','demo1')
demoFile=fullfile(getenv('MFILES'),'venus','mlechipsdoslDemos.mat');

% Other defaults
defval('stit','Chi-squared residuals')

% Normal use of the function
if ~ischar(Hk)
  
    % More default values
    defval('ah',krijetem(subnum(2,2)))
     
    % If "Hk" is a structure, assume it is a full "mleR" structure with
    % results from Matern parameter estimation for a data patch
    if isstruct(Hk)
      % If "thhat" exists, it is the figure's supertitle text
        if exist(thhat,'var')
            stit=thhat;
        end
        % Pull necessary values from the passed structure
        getfieldr(Hk,{'Hx' 'thhat' 'params' 'logli'});
    
        % Recover Fourier domain data and filter the wavenumbers
        Hk=reshape(tospec(Hx,params.NyNx),params.NyNx);
        k=knums(params);
        Hk(k>params.kiso)=NaN;
    end
    
    %% TAKE CARE OF A FEW THINGS UPFRONT
    
    % Degrees of freedom for 2X system (assuming a single field)
    df=2;
    % Bounds of X residuals to show (functions as coloraxis in panel 3)
    boundsX0=[0 3*df];
    % Bounds of 2X residuals to how (functions as axis bounds in panel 1)
    bounds2X0=2*boundsX0;
    % Color axis for the power spectral densities (panels 2 and 4)
    caxPSD=[-5 0];
    % Contours to be plotted on the power spectral density plots
    conPSD=[-5:-1];
    
    % Generate the 2D wavenumber axis for "Hk"
    [k,kx,ky]=knum2(params.NyNx,[(params.NyNx(1)-1)*params.dydx(1) ...
                     (params.NyNx(2)-1)*params.dydx(2)]);
        % Get order of magnitude of last wavenumber for scaling
    om=round(log10(ky(end)));
    
    % Set up the wavenumber/wavelength axes labels for panels 2-4
    xstr1=sprintf('x wavenumber (rad/m) %s10^{%i}','\times',om);
    ystr1=sprintf('y wavenumber (rad/m) %s10^{%i}','\times',om);
    xstr2='x wavelength (km)';
    ystr2='y wavelength (km)';
    
    % Prepare the overall figure
    fh=gcf;
    fig2print(fh,'landscape')
    set(fh,'PaperPositionMode','auto')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PANEL 1: TWICE CHI-SQUARED RESIDUALS AND QQ-PLOT [Top-left]
    axes(ah(1))
    
    % Get the spectral matrix used in evaluating "thhat"
    [~,~,~,~,~,Sb]=simulosl(thhat,params,0);
    
    % Calculate residuals, removing values at the zero wavenumber and
    % wavenumbers above the minimum wavenumber "params.kiso"
    Xk=abs(Hk(:)).^2./Sb(:);
    Xk(k==0)=NaN;
    Xk(k>params.kiso)=NaN;
    
    % Evaluate the loglihood (need scaled "thhat") and report it
    scl=10.^round(log10(abs(thhat)));
    L=logliosl(thhat./scl,params,Hk(:),k,scl);
    disp(sprintf('The loglihood is %8.3f',L))
    
    % Get binning info for the twice chi-squared residuals histogram
    binWidth=df/4;
    binCent=0.25*(1:2:49);
    [bdens,c]=hist(2*Xk(~isnan(Xk)),binCent);
    
    % Plot the histogram as a bar graph
    bdens=bdens/indeks(diff(c),1)/length(Xk(~isnan(Xk)));
    bb=bar(c,bdens,1);
    set(bb,'FaceC',grey)
    hold on
    
    % Plot the ideal chi-squared distribution
    refs=linspace(0,max(bounds2X0),100);
    plot(refs,chi2pdf(refs,df),'Linew',1.5,'Color','k')
    hold off
    
    % Labeling and cosmetic adjustments
    xlim(bounds2X0)
    xl1(1)=xlabel('quadratic residual 2X');
    ylim([0 max(bdens)*1.05])
    yl1(1)=ylabel('probability density');
    
    % Prepare a new axis for the quantile-quantile plot
    ah2(1)=laxis(ah(1),0,0);
    axes(ah2(1))
    
    % Obtain qq-plot data
    h=qqplot(2*Xk,ProbDistUnivParam('gamma',[df/2 2]));
    hx=get(h,'Xdata'); hx=hx{1};
    hy=get(h,'ydata'); hy=hy{1};
    delete(h)
    
    % Make the qq-plot
    qq0=plot(bounds2X0,bounds2X0,'k'); hold on
    qq=plot(hx,hy,'LineS','none','Marker','o','MarkerF','r',...
                    'MarkerE','r','MarkerS',2);
    
    % More cosmetic adjustments
    set(ah(1),'box','off')
    set(ah2(1),'xlim',get(ah(1),'xlim'),'yaxisl','r','box','off',...
               'xaxisl','t','color','none','ylim',get(ah(1),'xlim'))
    
    % Add labels and tickmarks
    ylr(1)=ylabel('quantile-quantile prediction');
    tkmks=bounds2X0(1):2*df:bounds2X0(2);
    set([ah(1) ah2(1)],'xtick',tkmks)
    set(ah2(1),'ytick',tkmks)
    
    % Add gridlines
    try
      gh=gridxy([4 8],[4 8],'LineStyle',':');
    end
    
    % For an info textbox, calculate the percent of 2X residuals
    % displayed, the maximum, mean, and variance of the 2X residuals
    tbstr{1}=sprintf('%5.2f%%',100*sum(binWidth*bdens(1:end-1)));
    tbstr{2}=sprintf('max(2X)=%5.2f',max(2*Xk));
    tbstr{3}=sprintf('mean(2X)=%5.2f',nanmean(2*Xk));
    tbstr{4}=sprintf('var(2X)=%5.2f',nanvar(2*Xk));
        
    % Calculate the mean of (X-df/2)^2 so it can be passed as output
    mag=nanmean((Xk-df/2).^2);
    tbstr{5}=sprintf('mean([X-%d]^2)=%5.2f',df/2,mag);
       
    % Give the x- and y-positions of the textbox
    tbx=repmat(bounds2X0(2)-bounds2X0(2)/30,1,length(tbstr));
    tby=[7.5 6.35-[0 1 2 3.5  ]];
    
    % Make the textbox(es) with unbordered fillboxes around them
    for i=1:length(tbstr)
        tb(i)=text(tbx(i),tby(i),tbstr{i},'HorizontalA','Right');
        gg=ext2lrtb(tb(i),[0.85],[0.9]); delete(tb(i)); hold on
        fb(i)=fillbox(gg+[0.8 0.8 0 0],'w');
        tb(i)=text(tbx(i),tby(i),tbstr{i},'HorizontalA','Right');
        hold off
        set(fb(i),'EdgeC','w','FaceC','w')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PANEL 2: PREDICTED 2D POWER SPECTRUM [Top-right]
    axes(ah(2))
    
    % The top-left and bottom-right wavenumbers, scaled
    c11=[kx(1) ky(end)]/10^om;
    cmn=[kx(end) ky(1)]/10^om;
    
    % Create a blurred wavenumber grid
    kb=knum2(params.blurs*params.NyNx,...
         [(params.blurs*params.NyNx(1)-1)*params.dydx(1) ...
         (params.blurs*params.NyNx(2)-1)*params.dydx(2)]);
    
    % Generate predicted, blurred power spectrum based on "thhat"
    Skpred=reshape(bluros(maternos(kb,thhat),params,1),params.NyNx);
    
    % Remove the zero wavenumber value (to avoid issues with the
    % colorscale) and those at wavenubmers above the isotropic cutoff
    Skpred(k==0)=NaN;
    Skpred(k>params.kiso)=NaN;
    
    % Scale to the largest predicted value and plot the power spectrum
    sclSkpred=log10(Skpred./max(Skpred(:)));
    imagefnan(c11,cmn,sclSkpred,'jet',caxPSD,[],[],0);
    
    % Label the wavenumber axes
    xl1(2)=xlabel(xstr1);
    yl1(2)=ylabel(ystr1);
    
    % For the wavelength axis, get tickmark values on the wavenumber axis
    xtk=get(ah(2),'xtick')*10^om;
    ytk=get(ah(2),'ytick')*10^om;
    
    % Convert to wavelengths (in m), recognizing the zero wavenumber
    xtkl=2*pi./xtk/1000;
        xtkl(isinf(xtkl))=[params.NyNx(2)-1]*params.dydx(2)/1000;
    ytkl=2*pi./ytk/1000;
        ytkl(isinf(ytkl))=[params.NyNx(1)-1]*params.dydx(1)/1000;
    
    % Create and label the wavelength axis
    [ah2(2),xl2(2),yl2(2)]=xtraxis(ah(2),xtk/10^om,round(xtkl),xstr2,...
        ytk/10^om,round(ytkl),ystr2);
    
    % Return to the main axis and prepare to plot contours
    axes(ah(2)); hold on
    
    % Calculate and plot contours
    [~,ch(1)]=contour(kx/10^om,fliplr(ky/10^om),sclSkpred,conPSD,'LineW',2);
    caxis(caxPSD)
    
    % Option to set the contours to black for visibility
    %set(hb,'color','k')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% PANEL 3: 2D X RESIDUALS [Bottom-left]
    axes(ah(3))
    
    % Plot the X residuals
    imagefnan(c11,cmn,reshape(Xk,params.NyNx),'gray',boundsX0,[],1);
    axis image
    
    if ~isnan(params.kiso)
      % Place a box giving the wavelength of the isotropic wavenumber cutoff
     [~,zy]=boxtex('ll',ah(3),sprintf('%2.0f km',2*pi./params.kiso/1000),...
                   12,[],[],1.1);
     set(zy,'fonts',get(zy,'fonts')-1)
    end
    
    % Label the wavenumber axes
    xl1(3)=xlabel(xstr1);
    yl1(3)=ylabel(ystr1);
    
    % Create and label the wavelength axis
    ah2(3)=xtraxis(ah(3),xtk/10^om,round(xtkl),xstr2,...
        ytk/10^om,round(ytkl),ystr2);
    
    % Remove the right y-label to avoid clutter
    delete(get(ah2(3),'ylabel'))
    
    axes(ah(3))
    
    % Add a colorbar
    [cb,xcb]=addcb('vert',boundsX0,boundsX0,'gray',df,1);
    set(xcb,'string','quadratic residual X')
    
    % Reposition the main axis
    set(ah(3),'position',[getpos(ah(3),1) getpos(ah(3),2) ...
        getpos(ah(2),3) getpos(ah(2),4)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PANEL 4: 2D TRUE POWER SPECTRUM [Bottom-right]
    axes(ah(4))
    
    % Convert (roughly) spectral domain topography to a power spectrum
    Sk=abs(Hk).^2;
    
    % Remove the zero wavenumber value (to avoid issues with the
    % colorscale) and those at wavenubmers above the isotropic cutoff
    Sk(k==0)=NaN;
    Sk(k>params.kiso)=NaN;
    
    % Scale to the largest predicted value and plot the power spectrum
    sclSk=log10(Sk./max(Skpred(:)));
    imagefnan(c11,cmn,sclSk,'jet',caxPSD,[],[],0);
    
    % Label the wavenumber axes
    xl1(4)=xlabel(xstr1);
    yl1(4)=ylabel(ystr1);
    
    % Create and label the wavelength axis
    ah2(4)=xtraxis(ah(4),xtk/10^om,round(xtkl),xstr2,...
        ytk/10^om,round(ytkl),ystr2);
    
    % Return to the main axis and prepare to plot contours
    axes(ah(4)); hold on
    
    % Place contours from the predicted power spectrum onto this true one
    [~,ch(2)]=contour(kx/10^om,fliplr(ky/10^om),sclSkpred,conPSD,'LineW',2);
    caxis(caxPSD)
    
    % Option to set the contours to black for visibility
    %set(hb,'color','k')
  
    
    %% FINAL COSMETIC ADJUSTMENTS
    
    % Adjust tickmarks
    longticks([ah ah2 cb])
    
    % Adjust the colorbar in panel 3
    axes(cb)
    set(cb,'YAxisLocation','right')
    moveh(cb,0.06)
    shrink(cb,1.3,1.27)
    
    % Adjust the main axes
    set(ah(2:4),'Box','On')
    
    % Give the overall figure a title
    axes(ah2(1))
    spt=title(ah2(1),stit);
        
    movev([ah ah2 cb],-.02)
    movev(spt,-.05)
    
    % Set figure background color to white
    set(gcf,'color','W','InvertH','off')
    
    % Collect output
    vars={mag,ah,ah2,cb,ch,spt};
    varargout=vars(1:nargout);

elseif strcmp(Hk,'demo1')
    % Set parameters for creation of a data patch
    fields={'dydx','NyNx','blurs'};
    defstruct('params',fields,{[20 20]*1e3,128*[1 1],2});
    % Random random parameters
    th0=max(round(rand(1,3).*[1 1 1]*10),[1 1 1])./[1e-4 1 1e-4];
    th0(2)=2+rand(1,1)*2;
    
    % Create the data patch, both in spatial and Fourier domain
    [Hx,~,params,k,Hk]=simulosl(th0,params); 
    
    % Set the isotropic wavenumber cutoff to some random fraction
    params.kiso=k(1,randi([params.NyNx(1)/2,params.NyNx(1)]));
    
    % Estimate the parameters via maximum-likelihood
    % Initial values are now handled well inside of it
    % thini=[10000 1 1e5];
    thini=[];
    [th,~,~,~,scl]=mleosl(Hx,thini,params);
    thhat=th.*scl;

    % Create a title
    stit=sprintf('est %s = [%i %5.2f %i]\ntru %s = [%i %5.2f %i]',...
                   '\theta',round(thhat(1)),thhat(2),round(thhat(3)),...
                   '\theta',round(th0(1)),th0(2),round(th0(3)));
     
    % Produce the residuals/power spectral density figure
    clf
    mlechipsdosl(reshape(Hk,params.NyNx),thhat,params,stit);
    
    % Plot the figure!
    figna=figdisp([],[],[],1);
    system(sprintf('epstopdf %s.eps',figna));
    system(sprintf('rm -rf %s.eps',figna));
    
elseif strcmp(Hk,'demo2')
  % Rename demo number/ID, and free up the input variables for reuse
    demoNum=Hk; Hk=[];
    defval('thhat',1)
    demoID=thhat; thhat=[];
    
    % If a string, "demoID" points toward a previously ran, saved demo
    if ischar(demoID)
        % Load the demo 2 structure from the demos savefile
        demo=modload(demoFile,demoNum);
        
        % If specific demo, load its variables
        if isfield(demo,demoID)
            fields={'Hk' 'params' 'thhat' 'kisos' 'pstit'};
            [Hk,params,thhat,kisos,pstit]=getfieldr(demo.(demoID),fields);
        
            % Generate the wavenumber grid
            k=knums(params);
            
        % If not, the requested demo does not exist
        else
            disp(sprintf('Requested %s %s cannot be found',demoNum,demoID))
        end
        
    % If not a string, do a fresh run of this demo
    else
      % For this condition, assume a custom field is being passed
        if exist('params','var')
            % Rename in the input variables to their proper variable names
            [Hx,params,pstit]=deal(params,stit,ah);
            % Convert the datafield to the spectral domain
            Hk=reshape(tospec(Hx,params.NyNx),params.NyNx);
            % Generate this variable for completeness when saving
            th0=[NaN NaN NaN];
            % Generate the wavenumber grid
            k=knums(params);
            % Otherwise, generate a random field
        else
            % Set parameters for creation of a data patch
            fields={'dydx','NyNx','blurs','kiso'};
            defstruct('params',fields,{[20 20]*1e3,128*[1 1],2,NaN});
            
            % Set the "true" data patch Matern parameter values
            th0=max(round(rand(1,3).*[1 1 1]*10),[1 1 1])./[1e-4 1 1e-4];
            th0(2)=2+rand(1,1)*2;
            
            % Create the data patch, both in spatial and Fourier domain
            [Hx,~,~,k,Hk]=simulosl(th0,params);
            Hk=reshape(Hk,params.NyNx);
            
            % Create a pretitle for late figures
            pstit=sprintf('%s_0=[%9.3e %6.3f %6.0f]\n %s_0=[%9.3e %6.3f %6.0f]',...
                   '\theta',th0(1),th0(2),th0(3));
        end
        
        keyboard
        
        % Get the isotropic cutoff wavenumbers at no corners and all
        % corners removed, and at approximate fourths in between
        kisos=k(1,round(params.NyNx(1)/2*(1+[1 0.75 0.5 0.25 0])));
        
        % Make a new "params" structure for varying the isotropic cutoff
        nparams=params;
        
        % For each isotropic cutoff, calculate the Matern parameters that
        % best describe the data field
        for i=1:length(kisos)
            nparams.kiso=kisos(i);
            [th,~,~,~,scl]=mleosl(Hx,[0.01 1 1e5],nparams);
            thhat(i,:)=th.*scl;
        end
        
        % For this condition, save the run as a new demo
        if demoID == 1
          % If the demos savefile doesn't exist, create it
            if exist(demoFile,'file') == 0
                readme='Stored data for MLECHIPSDOSL demos';
                save(demoFile,'readme');
            end
            
            % If a particular demo structure is in the demos savefile, load
            % it.  Otherwise, create an empty structure for the new demo
            if inMat(demoNum,demoFile)
                demo=modload(demoFile,demoNum);
            else
                demo=struct;
            end
            
            % Get identifying datestring for this demo, save the demo
            % results to the structure, and save the structure
            demo.(datestr(now,'mmmdd_yyyy'))=struct('Hx',Hx,'Hk',Hk,...
                'th0',th0,'params',params,'thhat',thhat,'kisos',kisos,...
                'pstit',pstit);
            modSave(demoFile,demoNum,demo,'-append')
        end
    end
    
    % Set a range of colors for plots
    clrs={'b' 'g' 'r' 'c' 'm'};

    % Find the x- and y-coordinates of the zero wavenumber
    [zx,zy]=find(~k);
    difer(k(zx,zy),[],[],NaN)
    
    % Get the 1D, radial wavenumber axis
    kk=k(zx,zx:end);
    
    % Get the "true" 1D, radially averaged power spectrum
    [~,Hk1d]=radavg(abs(Hk.^2));
    
    % Figure of how the 1D power spectrum changes as corners are cut off
    figure
    for i=1:length(kisos)
        
        % Plot the 1D power spectrum from the estimated Matern parameters
        semilogx(kk,maternos(kk,thhat(i,:)),clrs{i});
        hold on
        legtxt{i}=sprintf('kiso=%0.4g',kisos(i));
    end

    % Plot the "true" 1D power spectrum
    semilogx(kk,Hk1d,'ko-')
    hold on
    legtxt{end+1}='Rad. ave. spectrum';
    
    % Plot vertical lines showing where each wavenumber cutoff occurs
    for i=1:length(kisos)
        plot([kisos(i) kisos(i)],ylim,clrs{i})
    end
    
    % Annotate the figure
    legend(legtxt)
    xlabel('wavenumber')
    ylabel('power')
    title('1D power spectrum as corners are removed')
    
    % Figure of how the 1D power spectrum normalized by the estimated s2 
    % changes as corners are cut off
    figure
    for i=1:length(kisos)
        
        % Plot the 1D power spectrum from the estimated Matern parameters
        semilogx(kk,maternos(kk,thhat(i,:))/thhat(i,1),clrs{i});
        hold on
        legtxt{i}=sprintf('kiso=%0.4g',kisos(i));
    end
    
    % Plot vertical lines showing where each wavenumber cutoff occurs
    for i=1:length(kisos)
        plot([kisos(i) kisos(i)],ylim,clrs{i})
    end
    
    % Annotate the figure
    legend(legtxt{1:end-1})
    xlabel('wavenumber')
    ylabel('power')
    title('1D power spectrum/s2 as corners are removed')
    
    % Make series of MLECHIPSDOSL plots for each wavenumber cutoff
    for i=1:length(kisos)
        
        % Set the isotropic wavenumber cutoff
        params.kiso=kisos(i);
        
        % Create the supertitle for the figure
        stit=sprintf('%s\n %s=[%9.3e %6.3f %6.0f]; kiso=%0.4g',...
                       pstit,'\theta',thhat(i,1),thhat(i,2),thhat(i,3),kisos(i));
        
        % Make the figure
        figure  
        ah=krijetem(subnum(2,2));
        mlechipsdosl(Hk,thhat(i,:),params,stit,ah)
    end
end
