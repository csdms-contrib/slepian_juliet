function varargout=mlelcontosl(Hk,thhat,scl,params,covth,thRange,runinpar,stit,ah)
% [Lgrid,Lcon,thR,xcon,ycon,ah,spt]=MLELCONTOSL(Hk,thhat,scl,params,covth,thRange,runinpar,stit,ah)
%
% Calculates the log-likelihood surface in the parameter cross-plot space (3 2D
% parameter spaces, s2-nu, s2-rho, and nu-rho) for the provided data.
% The three calculated log-likelihood grids are centered on the estimate of the 
% two parameters described and hold the third parameter constant at its 
% estimated value; the grids are contoured a certain number of standard 
% deviations for each parameter. The standard deviation is calculated from the 
% covariance matrix provided. We expect the Matern estimate to be centered in a 
% log-likelihood valley.
%
% INPUT:
%
% Hk         The Fourier-domain data, e.g. from SIMULOSL
% thhat      The three estimated Matern parameters (scaled); optionally, provide
%            multiple thhats and NO Hk to calculate likelihood contours for 
%            many estimates 
% scl        The scaling of the Matern parameters and Hk (scl(1))
% params     Parameter set pertaining to the data, e.g. from SIMULOSL
% covth      Covariance matrix of the estimates, whichever way computed
% thRange    Boundary of desired figure; if specified, may reduce number of 
%            calculations required; [default extends -/+ 5 s.d. from thhat]
% runinpar     Flag for running each of the 2D parameter space loglihood
%            calculations in parallel [default]
% stit       Title for the outputted figure [defaulted]
% ah         A 4x4 group of axes handles [defaulted]
%
% OUTPUT
% 
% Lgrid      The grid of loglihood values; ordered to correspond with usual
%            MLEOSL pairings: s2-nu, s2-rh, nu-rh
% Lcon       The values to contour for, based on covth; order: s2-nu, s2-rh, nu-rh
% thR        Vectors defining the Matern parameter grids that the loglihoods
%            were evaluated for; order: s2-nu, s2-rh, nu-rh 
% xcon       Cell array for the horizontal axis values of each pairing
% ycon       Cell array for the vertical axis values of each pairing
% ah         Set of axis handles
% spt        Supertitle handle
% 
% 
% EXAMPLE:
%
% mlelcontosl('demo1')
%
% [Hx,th0,params,k,Hk]=simulosl;
% [thhat,covFHh,lpars,scl,thini,params,Hk,k]=mleosl(Hx,[],params,[],[],[],[],1,2);
% thRange=[thhat'.*scl'+2*sqrt(diag(covFHh{4})).*[-1 1]];
% thRange=[0.8*thhat' 1.2*thhat'];
% thRange = [];
% [Lgrid,Lcon,thR]=mlelcontosl(Hk,thhat,scl,params,covFHh{3},thRange,1);
%
% SEE ALSO:
%
% CONTOUR, LOGLIOSL
%
% Last modified by gleggers-at-alumni.princeton.edu, 05/26/2014
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2023
% Last modified by olwalbert-at-princeton.edu, 05/28/2025

% Default values for function inputs
defval('stit','Loglihood contours')
% Launch a figure that will not interfere with MLEOSL demo2 which calls this
% routine
figure(3); clf
defval('ah',krijetem(subnum(2,2)))
% Run in parallel?
defval('runinpar',1)
if runinpar
  NumWorkers=feature('numcores')-1;
  if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end
end

% Normal use of the function, i.e. not a demo
if ~ischar(Hk)
    % Default values for the loglihood grid and contours
    % Number of standard deviations the grid extends
    defval('stdOut',2);         
    % Standard deviations to which the contours refer
    stdCon=1:stdOut;
    % Fineness of the loglihood grids
    fine=50;

    % Default options for figure contruction
    % Shade the plot outside the first contour
    optShade=false;   
    % Label contour loglihood values
    optClabel=true;  
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
    
    % Scaling of the Matern parameters relative to the data vector
    scl2=scl; scl2(1)=1;
    % Scaling of the data vector
    scl3=[scl(1) 1 1];
    % The Matern estimate scaled relative to the data vector
    thhat2=thhat.*scl2;
    % The unscaled Matern estimate
    thhat3=thhat.*scl;

    % Find loglihood of the estimated Matern parameters
    logli=logliosl(K,thhat2,params,Hk);
    
    % Get the standard deviations of the estimated Matern parameters
    thhatsd=sqrt(diag(covth));
    % The standard deviations relative to the data vector
    thhatsd2=thhatsd./scl3';

    % Normalize the estimated Matern parameters covariance matrix
    covthn=covth./(sqrt(diag(covth))*sqrt(diag(covth)'));

    % Extract elements from the normalized covariance matrix
    % (cross-correlation coefficients)
    Csn=covthn(1,2);
    Csr=covthn(1,3);
    Cnr=covthn(2,3);
    
    % Set the Matern parameter range that the figures will span
    defval('thRange',[thhat2'+stdOut*thhatsd2.*[-1 1]]);
    
    % Calculate the Matern parameter values for each contour
    preCon=repmat(thhat2',1,length(stdCon))+thhatsd2*stdCon;
    
    % Set up vectors defining the Matern parameter grids on which
    % loglihoods will be evaluated    
    s2R=unique([linspace(thRange(1,1),thhat2(1),fine/2) ...
                linspace(thhat2(1),thRange(1,2),fine/2)]);
    nuR=unique([linspace(thRange(2,1),thhat2(2),fine/2) ...
                linspace(thhat2(2),thRange(2,2),fine/2)]);
    rhoR=unique([linspace(thRange(3,1),thhat2(3),fine/2) ...
                 linspace(thhat2(3),thRange(3,2),fine/2)]);
    
    % Combine all three vectors for easy reference
    thR=[s2R; nuR; rhoR];

    % This from MLEPLOS as opposed to Gabe's cycIndex
    pcomb=nchoosek(1:length(thhat),2);

    % Calculate loglihoods on Matern parameter grids. For each pairing of
    % parameters (s2-nu, s2-rho, and nu-rho - the third parameter being held
    % constant at its MLE estimated value), calculate the loglihood value,
    % stored in "preCon", that will be used for contouring the pairs of 
    % parameters on the loglihood grid defined by "thR"

    % Pre-allocate space for the loglihood grid
    Lgrid=zeros(fine-1,fine-1,3); 
    for i=1:length(thhat)
        % Find the pairwise combinations
        xi=pcomb(i,1); yi=pcomb(i,2);

      % Update on which loglihood grid is being constructed
      disp(sprintf('Constructing %s-%s loglihood grid, %s',...
                   thNam{xi},thNam{yi},datestr(now,15)));

      % Calculate each contour's loglihood with proposed th's scaled relative to
      % the data vector
      for j=1:length(preCon(xi,:))
        switch i
          case 1
            % s2-nu
            th=[preCon(xi,j) preCon(yi,j) thhat2(3)];
          case 2
            % s2-rho
            th=[preCon(xi,j) thhat2(2) preCon(yi,j)];
          case 3
            % nu-rho
            th=[thhat2(1) preCon(xi,j) preCon(yi,j)];
        end

         % Calculate and store contour loglihood value
         Lcon(i,j)=logliosl(K,th,params,Hk);
      end
      % Calculate loglihoods on the Matern parameter grid
      if runinpar==1
        for j=1:fine-1
          % j: horizontal axis 
          Lg=cell(1,fine-1);
          thRxij=thR(xi,j);
          switch i
            case 1
              parfor (kk=1:fine-1,NumWorkers)
                % kk: vertical axis 
                th=[thRxij thR(yi,kk) thhat2(3)];
                Lg{kk}=logliosl(K,th,params,Hk);
              end
              Lgrid(:,j,i)=cat(1,Lg{:});
            case 2
              parfor (kk=1:fine-1,NumWorkers)
                  th=[thRxij thhat(2) thR(yi,kk)];
                  Lg{kk}=logliosl(K,th,params,Hk);
              end
              Lgrid(:,j,i)=cat(1,Lg{:});
            case 3
              parfor (kk=1:fine-1,NumWorkers)
                th=[thhat(1) thRxij thR(yi,kk)];
                Lg{kk}=logliosl(K,th,params,Hk);
              end
              Lgrid(:,j,i)=cat(1,Lg{:});
          end
        end
      else
        % Not running in parallel; i: parameter pairing, j: horizontal, k:
        % vertical
        for j=1:length(thR(xi,:))
          for k=1:length(thR(yi,:))
            % Give a grid point's proposed theta, scaled relative to the data
            % vector
            switch i
              case 1
                th=[thR(xi,j) thR(yi,k) thhat2(3)];
              case 2
                th=[thR(xi,j) thhat2(2) thR(yi,k)];
              case 3
                th=[thhat2(1) thR(xi,j) thR(yi,k)];
            end
              
            % Calculate and store loglihood value for the grid point
            Lgrid(k,j,i)=logliosl(K,th,params,Hk);
          end
        end
      end
    end
    
    % Loglihood contour figure construction
    % Default values for figure aesthetics:  
    % Standard deviations to have tickmarks
    stds=-stdOut:stdOut; 

    % Get the text showing the scaling of each Matern parameter axis
    for i=1:length(thhat)
        scTxt{i}=sprintf('x 10^{%d}',round(log10(scl(i)))); 
    end

    % Calculate the values for tickmarks on the loglihood grid
    thTick=repmat(thhat3',size(stds))+thhatsd*stds;
    thTickLa=cell(size(thTick));

    % Convert those tickmark values to strings for tickmark labels
    for i=1:length(stds)
        % For tickmarks on even standard deviations, round and scale the value
        if mod(stds(i),2)==0
            for j=1:length(thhat)
                % Round to two decimal places
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
        % [xi,yi]=cycIndex(i,length(thhat));
        xi=pcomb(i,1); yi=pcomb(i,2);
        
        % Option to shade the plot background outside the largest contour
        if optShade
            imagefnan([thR(xi,1) thR(yi,1)],[thR(xi,end) thR(yi,end)],...
                      Lgrid(:,:,xi)-logli,[],[],[],[],...
                      max(abs(reshape(Lgrid(:,:,xi),1,[])-logli))/(Lcon(xi,1)-logli));
            hold on
        end
        
        % Place dotted gridlines at the center and "glp" standard deviations
        glp=stdOut;
        hold on
        glh(1)=plot([thhat3(xi) thhat3(xi)],...
                    [thR(yi,1) thR(yi,end)].*scl3(yi),'k:');
        glh(2)=plot([thhat3(xi) thhat3(xi)]+glp*thhatsd(xi),...
                    [thR(yi,1) thR(yi,end)].*scl3(yi),'k:');
        glh(3)=plot([thhat3(xi) thhat3(xi)]-glp*thhatsd(xi),...
                    [thR(yi,1) thR(yi,end)].*scl3(yi),'k:');
        
        glh(4)=plot([thR(xi,1) thR(xi,end)].*scl3(xi),...
                    [thhat3(yi) thhat3(yi)],'k:');
        glh(5)=plot([thR(xi,1) thR(xi,end)].*scl3(xi),...
                    [thhat3(yi) thhat3(yi)]+glp*thhatsd(yi),'k:');
        glh(6)=plot([thR(xi,1) thR(xi,end)].*scl3(xi),...
                    [thhat3(yi) thhat3(yi)]-glp*thhatsd(yi),'k:');

        % Plot the loglihood contours
        xcon{i}=thR(xi,:).*scl3(xi);
        ycon{i}=thR(yi,:).*scl3(yi);

        % FIGURE OUT how to replace negative parts of the grid with NaNs

        %Lgrid(xcon<=0,ycon<=0,xi)=NaN;
        [c,ch(i)]=contour(xcon{i},ycon{i},Lgrid(:,:,i),Lcon(i,:));
        
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
        thSol(i)=plot(thhat3(xi),thhat3(yi),'o');
        
        % Option for markers on the parameter pairs whose loglihoods gave the
        % plotted contours, i.e.thhat + 1, 2, 3 STDs
        if optSD
            hold on
            thCon(i)=plot(repmat(thhat3(xi),1,length(stdCon))+stdCon*thhatsd(xi),...
                          repmat(thhat3(yi),1,length(stdCon))+stdCon*thhatsd(yi),'o');
        end
        
        % Option for an "error ellipse" over the contour plot
        % Heretofore been unable to get this to work, but see MLEPLOS
        if optErrE == true
            covth2D=[covth(xi,xi) covth(xi,yi); covth(yi,xi) covth(yi,yi)];
            error_ellipse(covth2D,[thhat3(xi),thhat3(yi)]);
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
    strAnno{4}=sprintf('\\sigma^{2}=%6.4f %s %6.4f',thhat3(1),char(177),glp*thhatsd(1));
    strAnno{5}=sprintf('\\nu=%6.4f %s %6.4f',thhat3(2),char(177),glp*thhatsd(2));
    strAnno{6}=sprintf('\\rho=%6.4f %s %6.4f',thhat3(3),char(177),glp*thhatsd(3));
    strAnno{7}=sprintf('C_{\\sigma^{2}\\nu}=%5.4f',Csn);
    strAnno{8}=sprintf('C_{\\sigma^{2}\\rho}=%5.4f',Csr);
    strAnno{9}=sprintf('C_{\\nu\\rho}=%5.4f',Cnr);

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
    varns={Lgrid,Lcon,thR,xcon,ycon,ah,spt};
    varargout=varns(1:nargout);
elseif strcmp(Hk,'demo1')
    clear
    % Set parameters for creation of a data patch
    fields={'dydx','NyNx','blurs'};
    defstruct('params',fields,{[10 10]*1e3,64*[1 1],-1});
    % defstruct('params',fields,{[1 1]*1e0,64*[1 1],-1});
    % Random random parameters
    th0=max(round(rand(1,3).*[1 1 4]*10),[1 1 1])./[1e-4 1 1e-4];
    th0(2)=2+rand(1,1)*2;
    % Examples close to those we've already done
    th0=[1e6 2.5 2e4];
    % th0=[1 1.5 2];
    try
        % Create the data patch, both in spatial and Fourier domain
        [Hx,~,params,k,Hk]=simulosl(th0,params); 
        try
            % Estimate the parameters via maximum-likelihood
            [th,covFHhJ,lpars,scl,~,~,Hk]=mleosl(Hx,[],params,[],[],[],[],[],2);

            % Remind us where the loglihood was
            disp(sprintf('L = %6.2f',lpars{1}))

            % Produce the likelihood contours figure
            clf
            runinpar=1;
            [Lgrid,Lcon,thR]=mlelcontosl(Hk,th,scl,params,covFHhJ{4},[],runinpar);
            
            % Plot the figure! 
            figna=figdisp([],[],[],1);
        end
    end
end
