function varargout=mlelcontosl(Hk,thhat,params,covth,stit,ah)
% [lcont,ah,spt]=MLELCONTOSL(Hk,thhat,params,covth,stit,ah)
%
% For a Matern solution set (s2,nu,rho) describing a patch of data, this
% function shows in Matern parameter-space how the log-likelihood varies
% around the estimated solution. The function creates a set of three grids
% (s2-nu, nu-rho, and rho-s2) centered on the estimated set and extending a
% certain number of standard deviations in each parameter-direction. The
% loglihoods on that grid are calculated (the third parameter is held
% constant), and then that loglihood space is contoured using positive
% standard deviation pairs from the solution center. What should be seen is
% that the Matern solution set lies in a loglihood valley (as in finding the
% solution it was negative loglihood being minimized), and this can be used
% as a measure to see if the best Matern solution was found.
%
% INPUT:
%
% Hk         The Fourier-domain data, e.g. from SIMULOSL
% thhat      The three Matern parameters
% params     Parameter set pertaining to the data, e.g. from SIMULOSL
% covth      Covariance matrix of the estimates, whichever way computed
% stit       Title for the overall figure [defaulted]
% ah         A 4x4 group of axis handles [defaulted]
%
% OUTPUT
% 
% lcont      Likelihood contours
% ah         Set of axis handles
% spt        Supertitle handle
% 
% EXAMPLE:
%
% mlelcontosl('demo1')
% 
% SEE ALSO:
%
% CONTOUR, LOGLIOSL
%
% Last modified by gleggers-at-alumni.princeton.edu, 05/26/2014
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2023

% Default values for function inputs
defval('stit','Loglihood contours')
defval('ah',krijetem(subnum(2,2)))

% Normal use of the function, i.e. not a demo
if ~ischar(Hk)
    % Default values for the loglihood grid and contours
    % Number of standard deviations the grid extends
    stdOut=5;         
    % Standard deviations to which the contours refer
    stdCon=[1 2 3];
    % Fineness of the loglihood grids
    fine=50;

    % Default options for figure contruction
    % Shade the plot outside the first contour
    optShade=false;   
    % Label contour loglihood values
    optClabel=false;  
    % Show the STD pairs from which contours are calculated
    optSD=false;      
    % Plot an error ellipse about the calculated solution
    optErrE=false;
    
    % Names of the parameters, both simple and for text formatting
    thNam={'s2' 'nu' 'rho'};
    thForm={'\sigma^{2}' '\nu' '\rho'};

    % Loglihood grid calculation
    
    % Get wavenumber grid and number of wavenumbers used finding "thhat"
    K=knums(params);
    Kn=sum(~isnan(Hk));
    
    disp('FJS here watch that kiso is not counted we put in LOGLIOSL')
    
    % Get scaling of the estimated Matern parameters
    scl=10.^round(log10(abs(thhat)));

    % Find loglihood of the estimated Matern parameters
    logli=logliosl(K,thhat./scl,scl,params,Hk);
    
    % Get the standard deviations of the estimated Matern parameters
    thhatsd=sqrt(diag(covth));
    
    % Normalize the estimated Matern parameters covariance matrix
    covthn=covth./(sqrt(diag(covth))*sqrt(diag(covth)'));
    
    % Extract individual covariances from the normalized matrix
    Csn=covthn(1,2);
    Csr=covthn(1,3);
    Cnr=covthn(2,3);
    
    % Set the Matern parameter range the figures will encompass
    thRange=[thhat'-stdOut*thhatsd thhat'+stdOut*thhatsd];
    
    % Calculate the Matern parameter values for each contour
    preCon=repmat(thhat',1,length(stdCon))+thhatsd*stdCon;
    
    % Set up vectors defining the Matern parameter grids on which
    % loglihoods will be evaluated    
    s2R=unique([linspace(thRange(1,1),thhat(1),fine/2) ...
                linspace(thhat(1),thRange(1,2),fine/2)]);
    nuR=unique([linspace(thRange(2,1),thhat(2),fine/2) ...
                linspace(thhat(2),thRange(2,2),fine/2)]);
    rhoR=unique([linspace(thRange(3,1),thhat(3),fine/2) ...
                 linspace(thhat(3),thRange(3,2),fine/2)]);
    
    % Combine all three vectors for easy reference
    thR=[s2R; nuR; rhoR];
    
    % Calculate loglihoods on Matern parameter grids. For each pairing of
    % parameters (s2-nu, nu-rho, and rho-s2 - the third parameter being held
    % constant at its MLE estimated value), calculate the loglihood that
    % will be contours from the pairs of parameters in "preCon" as well as
    % the loglihood at points on the grid defined by "thR"
    for i=1:length(thhat)
        % Get the rotating indices
        [xi,yi]=cycIndex(i,3);
        
        % Update on which loglihood grid is being constructed
        disp(sprintf('Constructing %s-%s loglihood grid, %s',...
                     thNam{xi},thNam{yi},datestr(now,15)));
        
        % Calculate each contour's loglihood - this used to be in parallel
        for j=1:length(preCon(xi,:))
            % Give a loglihood contour's "estimated theta"
            switch xi
              case 1
                th=[preCon(xi,j) preCon(yi,j) thhat(3)];
              case 2
                th=[thhat(1) preCon(xi,j) preCon(yi,j)];
              case 3
                th=[preCon(yi,j) thhat(2) preCon(xi,j)];
            end
            
            % Scale the estimated theta
            scl=10.^round(log10(abs(th)));
            th=th./scl;

            % Calculate and store contour loglihood value
            Lcon(i,j)=logliosl(K,th,scl,params,Hk);
        end
        
        % Calculate loglihoods on the Matern parameter grid
        for j=1:length(thR(xi,:))
            for k=1:length(thR(yi,:))
                % Give a grid point's "estimated theta"
                switch xi
                  case 1
                    th=[thR(xi,j) thR(yi,k) thhat(3)];
                  case 2
                    th=[thhat(1) thR(xi,j) thR(yi,k)];
                  case 3
                    th=[thR(yi,k) thhat(2) thR(xi,j)];
                end
                
                % Scale the estimated theta
                scl=10.^round(log10(abs(th)));
                th=th./scl;
                
                % Calculate and store loglihood value for the grid point
                Lgrid(j,k,i)=logliosl(K,th,scl,params,Hk);
            end
        end
    end

    % Loglihood contour figure construction

    % Default values for figure aesthetics:  
    % Place to which to round axis tickmark labels
    rTo=-2;               
    % Standard deviations to have tickmarks
    stds=-stdOut:stdOut; 

    % Get the text showing the scaling of each Matern parameter axis
    for i=1:length(thhat)
        scTxt{i}=sprintf('x 10^{%d}',log10(scl(i))); 
    end

    % Calculate the values for tickmarks on the loglihood grid
    thTick=repmat(thhat',size(stds))+thhatsd*stds;
    thTickLa=cell(size(thTick));

    % Convert those tickmark values to strings for tickmark labels
    for i=1:length(stds)
        % For tickmarks on even standard deviations, round and scale the value
        if mod(stds(i),2)==0
            for j=1:length(thhat)
                % This used to be ROUND2
                % thTickLa{j,i}=sprintf('%0.2f',round2(thTick(j,i)/scl(j),rTo));
                thTickLa{j,i}=sprintf('%0.2f',round(thTick(j,i)/scl(j)*100)/100);
            end
            % For tickmarks on odd standard deviations, don't show any values
        else
            thTickLa(:,i)={''};
        end
    end

    % Special colors for when that used to work
    cols={'r' 'b' 'g'};

    % Make the loglihood grid contour plots on panels 1-3
    for i=1:length(thhat)
        % Activate the proper panel
        axes(ah(i))
        % Get the rotating indices
        [xi,yi]=cycIndex(i,length(thhat));
        
        % Option to shade the plot background outside the largest contour
        if optShade
            imagefnan([thR(xi,1) thR(yi,1)],[thR(xi,end) thR(yi,end)],...
                      Lgrid(:,:,xi)-logli,[],[],[],[],...
                      max(abs(reshape(Lgrid(:,:,xi),1,[])-logli))/(Lcon(xi,1)-logli));
            hold on
        end
        
        % Place dotted gridlines at the center and "glp" standard deviations
        glp=2;
        hold on
        glh(1)=plot([thhat(xi) thhat(xi)],[thR(yi,1) thR(yi,end)],'k:');
        glh(2)=plot([thhat(xi) thhat(xi)]+glp*thhatsd(xi),[thR(yi,1) thR(yi,end)],'k:');
        glh(3)=plot([thhat(xi) thhat(xi)]-glp*thhatsd(xi),[thR(yi,1) thR(yi,end)],'k:');
        
        glh(4)=plot([thR(xi,1) thR(xi,end)],[thhat(yi) thhat(yi)],'k:');
        glh(5)=plot([thR(xi,1) thR(xi,end)],[thhat(yi) thhat(yi)]+glp*thhatsd(yi),'k:');
        glh(6)=plot([thR(xi,1) thR(xi,end)],[thhat(yi) thhat(yi)]-glp*thhatsd(yi),'k:');

        % Plot the loglihood contours
        xcon=thR(xi,:);
        ycon=thR(yi,:);
        
        % FIGURE OUT how to replace negative parts of the grid with NaNs
        
        %Lgrid(xcon<=0,ycon<=0,xi)=NaN;
        % This is three contours but you're only getting to address all of them
        % Might need to go to individual contours
        [c,ch(i)]=contour(xcon,ycon,Lgrid(:,:,xi),Lcon(xi,:));
        
        % Option to label the contours with their loglihood values
        if optClabel
            clabel(c,ch(i))
        end
        
        % Get the handles to the contours and reverse their order - the contours
        % are naturally in ascending order of negative loglihood (or, rather,
        % decreasing order of loglihood) - this used to be for a different version
        % conts=getkids(ch(i));
        % conts=conts(end:-1:1);
        % Get the values of the contours of this panel
        % [~,lcont]=getcontour(c); 
        %lcont=lcont(:);
        % Grab the unique values (should only be three) and set colors for them
        % ucont=unique(lcont);
        % Check that these assumptions are satisfied
        % difer(length(ucont)-length(cols),[],[],NaN)
        % difer(length(lcont)-length(conts),[],[],NaN)
        % Set a color for each unique contour - for an older version
        %for j=1:length(lcont)
        %    set(conts(j),'EdgeColor',cols{ucont==lcont(j)})
        %end
        
        % Plot a center point showing the estimated Matern solution
        hold on
        thSol(i)=plot(thhat(xi),thhat(yi),'o');
        
        % Option for markers on the parameter pairs whose loglihoods gave the
        % plotted contours, i.e.thhat + 1, 2, 3 STDs
        if optSD
            hold on
            thCon(i)=plot(repmat(thhat(xi),1,length(stdCon))+stdCon*thhatsd(xi),...
                          repmat(thhat(yi),1,length(stdCon))+stdCon*thhatsd(yi),'o');
        end
        
        % Option for an "error ellipse" over the contour plot
        % Heretofore been unable to get this to work, but see MLEPLOS
        if optErrE == true
            covth2D=[covth(xi,xi) covth(xi,yi); covth(yi,xi) covth(yi,yi)];
            error_ellipse(covth2D,[thhat(xi),thhat(yi)]);
        end
        
        % Adjust the axis
        axis square tight
        
        % Set the x- and y-axis ticks and tick labels
        set(ah(i),'XTick',thTick(xi,:))
        set(ah(i),'XTickLabel',thTickLa(xi,:))
        
        set(ah(i),'YTick',thTick(yi,:))
        set(ah(i),'YTickLabel',thTickLa(yi,:))
        
        % Show axes scales (provided that scale is not 10^0)
        if scl(xi) ~= 1
            ttx(i)=text(thTick(xi,end-1),...
                        thTick(yi,1)-1.5*thhatsd(yi),scTxt(xi));
            set(ttx(i),'FontSize',max(get(gca,'FontSize')-2,6))
        end
        
        if scl(yi) ~= 1
            tty(i)=text(thTick(xi,1)-1.5*thhatsd(xi),...
                        thTick(yi,end-1)+1.5*thhatsd(yi),scTxt(yi));
            set(tty(i),'FontSize',max(get(gca,'FontSize')-2,6))
        end
        
        % Set the x- and y-axis labels and title
        xl(i)=xlabel(sprintf('%s',thForm{xi}));
        yl(i)=ylabel(sprintf('%s',thForm{yi}));
        
        tt=title(sprintf('log-likelihood in (%s,%s)',thForm{xi},thForm{yi})); 
    end

    % Switch to the fourth panel of the image for annotations
    axes(ah(4))

    % Set the text of the annotations
    strAnno{1}=sprintf('K=%5i',Kn);
    strAnno{2}=sprintf('k_{iso}=%0.4g',params.kiso);
    strAnno{3}=sprintf('L=%7.3f',-logli);
    strAnno{4}=sprintf('\\sigma^{2}=%6.4f %s %6.4f',thhat(1),char(177),glp*thhatsd(1));
    strAnno{5}=sprintf('\\nu=%6.2f %s %6.2f',thhat(2),char(177),glp*thhatsd(2));
    strAnno{6}=sprintf('\\rho=%5.0f %s %5.0f',thhat(3),char(177),glp*thhatsd(3));
    strAnno{7}=sprintf('C_{\\sigma^{2}\\nu}=%5.2f',Csn);
    strAnno{8}=sprintf('C_{\\sigma^{2}\\rho}=%5.2f',Csr);
    strAnno{9}=sprintf('C_{\\nu\\rho}=%5.2f',Cnr);

    % Get x- and y-positions for annotations
    xpos=get(ah(4),'XLim');
    ypos=get(ah(4),'YLim');

    % Set annotations
    for i=1:length(strAnno)
        ant(i)=text(xpos(2)*0.3,ypos(2)*0.9-(i-1)*ypos(2)/8,strAnno{i});
    end

    % Turn off the axis
    axis off

    % Various cosmetic adjustments
    fig2print(gcf,'portrait')
    longticks(ah(1:3))
    set(ah(1:3),'xgrid','off','ygrid','off','box','on')
    set(thSol(:),'MarkerF','k','MarkerE','r','MarkerS',4)
    set(ch(:),'LineW',1)
    serre(ah(1:2),[],'across')
    serre(ah(3:4),[],'across')

    % Make a supertitle for the whole figure
    spt=supertit(ah,stit);
    movev(spt,0.48)

    % Set figure background color to white
    set(gcf,'color','W','InvertH','off')

    % Collect output
    vars={lcont,ah,spt};
    varargout=vars(1:nargout);
elseif strcmp(Hk,'demo1')
    clear
    % Set parameters for creation of a data patch
    fields={'dydx','NyNx','blurs'};
    defstruct('params',fields,{[10 10]*1e3,128*[1 1],-1});
    % Random random parameters
    th0=max(round(rand(1,3).*[1 1 4]*10),[1 1 1])./[1e-4 1 1e-4];
    th0(2)=2+rand(1,1)*2;
    % Examples close tho those we've already done
    th0=[1e6 2.5 2e4];

    try
        % Create the data patch, both in spatial and Fourier domain
        [Hx,~,params,k,Hk]=simulosl(th0,params); 
        
        % imagesc(v2s(Hx,params)); axis image

        try
            % Estimate the parameters via maximum-likelihood
            thini=[];
            [th,covFHh,lpars,scl]=mleosl(Hx,thini,params);
        end
        % Scale up the outcome
        thhat=th.*scl;
        
        % Remind us where the loglihood was
        disp(sprintf('L = %6.2f',lpars{1}))
        
        % Produce the likelihood contours figure
        clf
        mlelcontosl(Hk,thhat,params,covFHh{3});
        
        % Plot the figure! 
        figna=figdisp([],[],[],2);
    end
end
