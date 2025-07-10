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
% thRange    Boundary of desired figure in parameter space [defaulted]
% runinpar   Flag for running each of the 2D parameter space loglihood
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
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2025
% Last modified by olwalbert-at-princeton.edu, 06/23/2025

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
    stdCon=1:1:stdOut;
    stdCon=0.2:0.2:stdOut
    % Fineness of the loglihood grids
    fine=200; fine=50

    % Default options for figure contruction
    % Shade the plot outside the first contour
    optShade=false;   
    % Label contour loglihood values
    optClabel=true;  
    % Show the STD pairs from which contours are calculated
    optSD=true;      
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
    defval('thRange',[repmat(thhat2',1,2)+stdOut*thhatsd2*[-1 1]]);
    
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

    keyboard
    
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
elseif strcmp(Hk,'demo2')
  % Calculate a 3D loglikelihood ellipsoid and volume for a dataset very close 
  % to th0 and compare to empirical and predicted ellipsoids and volumes

  % Load an experiment
  datum='27-May-2025'; trims=100;
  [th0,thhats,p,covX,covavhs,thpix,~,~,~,~,momx,covXpix,covF0]=osload(datum,trims);
  mobs=mean(thhats);
  cobs=cov(thhats);
  ncobs=cobs./[diag(sqrt(cobs))*diag(sqrt(cobs))']';
  sobs=std(thhats);

  n='mleosl';
  f4=sprintf('%s_diagn_%s',n,datum);
  np=3;
  % The loglikelihoods from estimate
  [~,~,~,~,L]=osrdiag(f4,pwd,np); 

  % Predict the covariance at th0
  p.taper=1;
  covth=covthosl(th0,p,2,[1 1 1]);

  % Go fish... for a data vector whose estimate is within P% of th0
  P=0.02;
  thhat=NaN;scl=NaN;
  while any(abs(thhat.*scl-th0)./th0>P) | isnan(thhat)
     % Simulate from the mean observation
     ps=p; ps.blurs=Inf; Hx=simulosl(th0,ps); 

     % This is where it would be good to load the simulation parameters; in
     % case we came in from a simulation with p.blurs=-1, this will need to
     % be uncommented:
     % p.blurs=-1; Hx=simulosl(mobs,p);

     % Calculate the MLE for the data vector
     [thhat,~,lpars,scl,~,~,Hk]=mleosl(Hx,[],p);
  end
  % Calculate the ellipsoids for the observed and predicted covariance
  e{1}=hyperellipsoid(mobs,cobs,3,[]);
  [e{2},V]=hyperellipsoid(th0,covth,3,[]);
  % Calculate the loglikelihoods for the parameter space (volume) straight from
  % LOGLIOSL
  thR=[minmax(thhats(:,1));minmax(thhats(:,2));minmax(thhats(:,3))].*[0.95 1.05];
  fine=40;
  thRlin{1}=linspace(thR(1,1),thR(1,2),fine);
  thRlin{2}=linspace(thR(2,1),thR(2,2),fine);
  thRlin{3}=linspace(thR(3,1),thR(3,2),fine);
  [thg{1},thg{2},thg{3}]=meshgrid(thRlin{1},thRlin{2},thRlin{3});
  thg{1}=thg{1}(:);
  thg{2}=thg{2}(:);
  thg{3}=thg{3}(:);
  K=knums(p);
  parfor kk=1:fine^3
    Lg{kk}=logliosl(K,[thg{1}(kk)./scl(1) thg{2}(kk) thg{3}(kk)],p,Hk);
  end
  thg{1}=reshape(thg{1},[fine fine fine]);
  thg{2}=reshape(thg{2},[fine fine fine]);
  thg{3}=reshape(thg{3},[fine fine fine]);
  Lgrid=reshape(cat(1,Lg{:}),[fine fine fine]);

  % 95% ellipsoid for the observed covariance colormapped by the loglikelihood
  % calculated along its surface for a data vector that is within P% of the truth
  es{1}=e{1}(:,:,1);es{1}=es{1}(:);
  es{2}=e{1}(:,:,2);es{2}=es{2}(:);
  es{3}=e{1}(:,:,3);es{3}=es{3}(:);
  parfor kk=1:size(es{1},1)
    Le{kk}=logliosl(K,[es{1}(kk)./scl(1) es{2}(kk) es{3}(kk)],p,Hk);
  end
  Lc{1}=reshape(cat(1,Le{:}),size(e{2},[1 2]));

  % 95% predicted ellipsoid colored by loglikelihood
  es{1}=e{2}(:,:,1);es{1}=es{1}(:);
  es{2}=e{2}(:,:,2);es{2}=es{2}(:);
  es{3}=e{2}(:,:,3);es{3}=es{3}(:);
  parfor kk=1:size(es{1},1)
    Le{kk}=logliosl(K,[es{1}(kk)./scl(1) es{2}(kk) es{3}(kk)],p,Hk);
  end
  Lc{2}=reshape(cat(1,Le{:}),size(e{2},[1 2]));

  % Create a figure with ellipsoids in the top row: emp, pred, like
  % and density heat maps in the bottom row: emp, pred, like
  [ah,ha,H]=krijetem(subnum(2,3));
  % Discretize the colormap so that we can see some contouring of the surfaces
  colormap(pink(18))
  % FaceAlpha for the surfaces
  fa=0.9;
  axes(ah(1))
  % 95% ellipsoid for thhats -- parameterized by mean(thhats) and
  % cov(thhats)
  s{1}=surf(e{1}(:,:,1),e{1}(:,:,2),e{1}(:,:,3),'FaceColor','k','FaceAlpha',0.05,...
         'EdgeColor','none');
  hold on
  s{2}=surf(e{2}(:,:,1),e{2}(:,:,2),e{2}(:,:,3),'FaceColor','b','FaceAlpha',0.05,...
         'EdgeColor','none');
  h(1)=scatter3(thhats(:,1),thhats(:,2),thhats(:,3),10,L,...
               'filled','MarkerFaceAlpha',fa);
 % [X,Y,Z]=ellipsoid(mobs(:,1),mobs(:,2),mobs(:,3),...
 %                   numstd*sobs(:,1),numstd*sobs(:,2),numstd*sobs(:,3));
 % s=surf(X,Y,Z,'FaceAlpha',0.2);
 % e=eig(ncobs);
 % rotate(s,[0,0,1],e(1));
 % rotate(s,[0,1,0],e(2));
 % rotate(s,[0,0,1],e(3));
 % rotate(s,[-1,0,0],atand(ncobs(2,3)));
 % rotate(s,[0,-1,0],atand(ncobs(1,3)));
 % rotate(s,[0,0,-1],atand(ncobs(1,2)));
  % p(1)=plot3();
  %lg=legend(sprintf('N=%i',size(thhats,1)))
  ti(1)=title({'\bf{Empirical}','estimates'});
  hold off

  axes(ah(2))
  s{3}=surf(e{2}(:,:,1),e{2}(:,:,2),e{2}(:,:,3),Lc{1},...
           'EdgeColor','none','FaceAlpha',fa);
%  s{3}=surf(e{2}(:,:,1),e{2}(:,:,2),e{2}(:,:,3),'FaceColor','k','FaceAlpha',0.1,...
%         'EdgeColor',clrs(2,:),'EdgeAlpha',0.05);
  ti(2)=title({'\bf{Empirical}','95% ellipsoid'});

  axes(ah(3))
  s{4}=surf(e{2}(:,:,1),e{2}(:,:,2),e{2}(:,:,3),Lc{2},...
           'EdgeColor','none','FaceAlpha',fa);
  ti(3)=title({'\bf{Predicted}','95% ellipsoid'});

  axes(ah(4))
  % Cross-section cut down rho revealing interior in s2 and nu planes; easiest
  % way is to set chunk to be cut to NaN and set the AlphaData for the plot
  Lgrid4=Lgrid;
  Lgrid4(thg{1}<th0(1)&thg{2}<th0(2))=deal(NaN);
  hs{1}=slice(thRlin{1},thRlin{2},thRlin{3},Lgrid4,...
              [minmax(thRlin{1}) min(thg{1}(thg{1}>th0(1)))],...
              [minmax(thRlin{2}) min(thg{2}(thg{2}>th0(2)))],...
              minmax(thRlin{3}));
  set(hs{1},'EdgeColor','none','FaceAlpha',fa,'FaceColor','interp');

  % th0slc=[th0(1) max(thRlin{1});...
  %         th0(2) max(thRlin{2});...
  %         th0(3) max(thRlin{3})];
  % slice(thRlin{1},thRlin{2},thRlin{3},Lgrid,...
  %       th0slc(1,1),th0slc(2,1),th0slc(3,1),'FaceAlpha',0.2)
  % hold on
  % slice(thRlin{1},thRlin{2},thRlin{3},Lgrid,th0slc(1,2),th0slc(2,2),th0slc(3,2))

  axes(ah(5))
  % Cross-section cut through along nu revealing interior in s2 and rh planes
  Lgrid5=Lgrid;
  Lgrid5(thg{1}<th0(1)&thg{3}>th0(3))=deal(NaN);
  hs{2}=slice(thRlin{1},thRlin{2},thRlin{3},Lgrid5,...
              [minmax(thRlin{1}) min(thg{1}(thg{1}>th0(1)))],...
              minmax(thRlin{2}),...
              [minmax(thRlin{3}) max(thg{3}(thg{3}<th0(3)))]);
  set(hs{2},'EdgeColor','none','FaceAlpha',fa,'FaceColor','interp');

  cb=colorbar('horizontal');cb.Label.String='loglikelihood';
  cb.Location='southoutside';
  movev(cb,-0.06)

  % Make sure that adding the colorbar did not permanently shrink the image axis
  ah(5).Position(4)=ah(4).Position(4);
  ah(5).Position(2)=ah(4).Position(2);

  axes(ah(6))
  % Planar cut along the eigenvector of the long axis of the error ellipsoid
  s2s=e{2}(:,:,1);s2s=s2s(:);nus=e{2}(:,:,2);nus=nus(:);rhs=e{2}(:,:,3);rhs=rhs(:);
  % rotang=rad2deg(1/2*atan2(2*covth(2,3),covth(2,2)-covth(3,3)));
  % s2-rh, rotated about the eigenvector for the nu axis
  rotang{1}=rad2deg(atan2(unique(rhs(s2s==max(s2s))),max(s2s)));
  hslice=surf(thRlin{1},thRlin{2},zeros(fine));
  %rotate(hslice,-1.*V(:,3),rotang{1})
  % or just rotate about nu to prevent jagged edges
  rotate(hslice,[0 -1 0],rotang{1})
  xd{1}=get(hslice,'XData');
  yd{1}=get(hslice,'YData');
  zd{1}=get(hslice,'ZData');
  delete(hslice)
  % nu-rh, rotated about the s2 axis
  % rotang{3}=rad2deg(atan2(unique(rhs(nus==max(nus))),max(nus)));
  % hslice=surf(thRlin{1},thRlin{2},zeros(fine));
  % rotate(hslice,[1 0 0],rotang{3})

  hs{3}=slice(thRlin{1},thRlin{2},thRlin{3},Lgrid,xd{1},yd{1},th0(3)+zd{1});
  hs{3}.FaceColor='interp';
  hs{3}.EdgeColor='none';
  hs{3}.DiffuseStrength=0.8;
  
  axes(ah(6))
  hold on 
  hx=slice(thRlin{1},thRlin{2},thRlin{3},Lgrid,max(thRlin{1}),[],[]);
  hx.FaceColor='interp';
  hx.FaceAlpha=fa;
  hx.EdgeColor='none';
  hy=slice(thRlin{1},thRlin{2},thRlin{3},Lgrid,[],minmax(thRlin{2}),[]);
  set(hy(1),'FaceColor','interp','FaceAlpha',fa/3,'EdgeColor','none')
  set(hy(2),'FaceColor','interp','FaceAlpha',fa/3,'EdgeColor','none')
  hz=slice(thRlin{1},thRlin{2},thRlin{3},Lgrid,[],[],min(thRlin{3}));
  hz.FaceColor='interp';
  hz.FaceAlpha=fa;
  hz.EdgeColor='none';

  for ind=1:6
    axes(ah(ind))
    % plot th0 
    hold on
    scatter3(th0(:,1),th0(:,2),th0(:,3),'ko','filled')
    if ind>1
      % plot the estimate associated with the data vector we use to calculate
      % the loglikelihoods
      scatter3(thhat(:,1).*scl(:,1),thhat(:,2).*scl(:,2),thhat(:,3).*scl(:,3),'c^','filled')
    end
    % set x,y,z,c limits to be the same
    xlim(minmax([minmax(thRlin{1}) minmax(e{1}(:,:,1)) minmax(e{2}(:,:,1))]))
    ylim(minmax([minmax(thRlin{2}) minmax(e{1}(:,:,2)) minmax(e{2}(:,:,2))]))
    zlim(minmax([minmax(thRlin{3}) minmax(e{1}(:,:,3)) minmax(e{2}(:,:,3))]))
    caxis(minmax(L))
    % other cosmetics
    xlabel('\sigma^2');  ylabel('\nu');  zlabel('\rho')
    grid on
    box on
    longticks
    % The azimuth and elevation perspective
    azel=[-50 7];
    view(azel)
  end
  movev([ah(1:3)],-0.04)
  movev([ah(4:6)],0.04)

  % Some additional 3D visualization ideas for another time:
  %
  % e{3}=isosurface(thg{1},thg{2},thg{3},Lgrid,12);
  % s(3)=patch(e{3});
  %
  % scatter3(thg{1}(:),thg{2}(:),thg{3}(:),10,Lgrid(:))
  %
  % surf(thRlin{1},thRlin{2},Lgrid(:,:,5),Lgrid(:,:,5))
  %
  % contourslice(Lgrid,[1:2:fine],[],[],10)

  % Information about the experiment
  sti=sgtitle(sprintf('th0=[%0.2e %0.2f %0.2e]; %ix%i grid; %ix%i km; blur %s; %s',...
        th0,p.NyNx,p.NyNx.*p.dydx./1000,ps.blurs));
  linkprop(ah,'CameraPosition')
  keyboard
  set(gcf,'Renderer','painters')
  figna=figdisp([],sprintf('demo2_%s',date),[],1);
  saveas(gcf,sprintf('mlelcontosl_demo2_%s.fig',date),'fig')
end
