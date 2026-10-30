function samplingcases(th,params,est,xver,fixflag)
% SAMPLINGCASES(TH,P,EST,XVER,FIXFLAG)
%
% How does sampling pattern impact Matern parameter estimation?
%
% Input:
%
% th        Matern parameters
% params    Grid parameters
% est       Select a figure style
%           0 - Create a figure in the style of Fig 1, ie, only simulated fields
%          [1 - Create a figure in the style of Fig 3, ie, generate parameter 
%               estimates with variance calculations for 100 realizations;
%               estimates will be generated and stored with a unique hash a la 
%               MLEOSL('DEMO2')]-- possibly to be removed
%           2 - Create a figure in the style of Fig 3 NEW
%           3 - Create a figure comparing a full field, a padded rectangle, and
%           a smaller rectangle
% xver      Number of iterations, borrowing naming convention from BLUROSY
% fixflag   Flag to set the fixed case to the WRONG nu, 5/2; otherwise fix to
%           the correct value of nu that we know from our simulation
%
% Examples:
%
% samplingcases([],[],0) 
% samplingcases([],[],2,[]) 
% samplingcases([],[],3,[])
%
% Output:
%
% <>        Generates and saves a figure (and possibly runs and stores many estimate
%           calculations)
%
% Last modified by olwalbert-at-princeton.edu 07/06/26

% Flag for Figure 1 or Figure 3
defval('est',0)
% Field size 
defval('params',[188 233]+randi(20,[1 2]))
% Rely on default scaling defined within MUCKIT or MASKIT
defval('scl',[]);
% Some combinations - SIMULOSL
%defval('th',[1 1.5 5e3]);
defval('th',[1 1.5 5]);
p.NyNx=[params(1) params(2)];
p.dydx=[1 1];
%p.dydx=1e3*[1 1];
p.blurs=Inf;
defval('fixflag',0)

datum=date;
%datum='30-Aug-2024';
%p.NyNx=[153 178]; % est==0 dims for figure 1
%p.NyNx=[188 253]; % est==3 dims for figure 3
%datum='25-Jun-2024'; %france mask for est==2; for original est==3;
%datum='21-Jun-2024'; %duck mask
%datum='20-Jun-2024'; 

if est==0
    % Create a figure in the style of Fig 1 of the sampling and special cases
    % paper, demonstrating the sampling cases as simulations only; many of the 
    % defaults and inputs of interest mirror DEMO2
    
    % We will show a panel for a (1) full rectangular field, (2) a field within
    % a geographic boundary, (3) speckled field, (4) an interior rectangular
    % subfield, (5) a field exterior to the geographic boundary ('anti-mask'),
    % and (6) the complement of the speckled field ('anti-muck'). We will count
    % down the column by defining the axes by the 'ha' index
    figure(1)
    clf
    [ah,ha,H]=krijetem(subnum(2,3));
    %Hx=simulosl(th,p);
    p.NyNx=[161 180]; % est==0 dims for figure 1
    for hnd=1:3
        if hnd==1
            % The internal and external rectangular case
            Hx=simulosl(th,p);
            p.taper=zeros(p.NyNx);
            % for only the inner 25% of the grid:
            qy=floor(p.NyNx(1)/4);qx=floor(p.NyNx(2)/4);
            cy=floor(p.NyNx(1)/2);cx=floor(p.NyNx(2)/2);
            p.taper(cy-qy:cy+qy,cx-qx:cx+qx)=deal(1.);
            % for the inner 75% of the grid:
            %qy=floor(p.NyNx(1)/16);qx=floor(p.NyNx(2)/16);
            %cy=floor(p.NyNx(1)/2);cx=floor(p.NyNx(2)/2);
            %p.taper(cy-7*qy:cy+7*qy,cx-7*qx:cx+7*qx)=deal(1.);
            p.mask=p.taper; 
            [v,I,w,vw,~,cr]=maskit(Hx,p,scl,Hx);
            % sum(I(:))/prod(p.NyNx)
            top=v;
            p2=p;p2.taper=~p.taper;p2.mask=~p.mask;
            [v,I,w,vw,~,cr]=maskit(Hx,p2,scl,Hx);
            bottom=v;mt=0;
        elseif hnd==2
            % The geographic boundary case; 
            %p.mask='france';[~,~,I]=maskit(rand(p.NyNx),p);p.taper=I;
            %p.mask='duck';[~,~,I]=maskit(rand(p.NyNx),p);p.taper=I;
            p.taper=1;
            p.mask='france';
            Hx=simulosl(th,p);
            scl=0.44; % this gives a mask retaining about 25% of the total grid
            % scl=1.42; % this gives a mask retaining about 75% of the grid --
            % this will obviously cut off the boundaries...
            [v,I,w,vw,scl,cr]=maskit(Hx,p,scl,Hx);
            % sum(I(:))/prod(p.NyNx)
            top=w;bottom=Hx;bottom(~~I(:))=NaN;mt=1;
        elseif hnd==3
            % The speckle case
            scl=0.25;
            %scl=0.75;
            p.taper=1;
            p.mask='random';
            Hx=simulosl(th,p);
            [v,I,w,scl,cr]=muckit(Hx,p,scl);
            %[~,I]=muckit(rand(p.NyNx),p,scl);p.taper=~I;
            top=v;bottom=w;mt=0;
        else                                                                    
            error('there is an indexing error; debug creation of subfigures')   
        end                                                                     
        axes(ha((hnd-1)*2+1))                                                   
        [h,t]=plotit(top,p,cr,th,mt);
        longticks
        set(h.Parent,'XTicklabel',[]);
        if hnd>1
            set(h.Parent,'YTickLabel',[]);
        end
        axes(ha((hnd-1)*2+2))                                                   
        [h2,t2]=plotit(bottom,p,cr,th,0);
        longticks
        if hnd>1
            set(h2.Parent,'YTickLabel',[]);
        end
    end  
    ho2=getpos(ha);o2=serre(H',1,'down');
    o3=serre(H,1/2,'across');
    movev(H,0.05);
    movev(H(2,:),0.035);
    figdisp([],sprintf('samplingfig1_%3.3i_%3.3i_fld',p.NyNx),[],1);
elseif est==1
    % Number of iterations
    defval('xver',100);
    % Method of computation
    defval('method','efs')

    % Once we find a grid-size we like through experimentation, let's set this
    % size for future program calls
    p.NyNx = [188, 233]; th = [1 1.5 5];

    % SET THE CORRECT P.TAPER
    filesuff=append('rect1_',datum);
    p_r1=p;
    try
        thhat_r1=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r1,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_r1,[],filesuff);
        thhat_r1=load(append('mleosl_thhat_',filesuff));
        while length(thhat_r1)<xver 
            mleosl('demo1',xver,th,p_r1,[],filesuff);
            thhat_r1=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    try
        thhat_r1_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_r1_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_r1,[],filesuff,th,[1 0 1]);
        thhat_r1_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r1_f,1)<xver
            mleosl('demo1',xver,th,p_r1,[],filesuff,th,[1 0 1]);
            thhat_r1_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The outer-rectangle case
    p_r2=p;
    p_r2.taper=ones(p_r2.NyNx);
    qy=floor(p_r2.NyNx(1)/4);qx=floor(p_r2.NyNx(2)/4);
    p_r2.taper(qy:3*qy,qx:3*qx)=deal(0.);
    p_r2.mask=p_r2.taper; 
    filesuff=append('rect2_',datum);
    try
        thhat_r2=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r2,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_r2,[],filesuff);
        thhat_r2=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r2,1)<xver
            mleosl('demo1',xver,th,p_r2,[],filesuff);
            thhat_r2=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    try
        thhat_r2_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_r2_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_r2,[],filesuff,th,[1 0 1]);
        thhat_r2_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r2_f,1)<xver 
            mleosl('demo1',xver,th,p_r2,[],filesuff,th,[1 0 1]);
            thhat_r2_f=load(append('mleosl_thhat_',filesuff));
        end
    end


    % The mask interior case
    p_m1=p;p_m1.mask='duck';
    %p_m1=p;p_m1.mask='ukraine';
    [~,I]=maskit(ones(p_m1.NyNx),p_m1);
    p_m1.taper=I;
    filesuff=append('mask1_',datum);
    try
        thhat_m1=load(append('mleosl_thhat_',filesuff));
        if size(thhat_m1,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_m1,[],filesuff);
        thhat_m1=load(append('mleosl_thhat_',filesuff));
        while size(thhat_m1,1)<xver 
            mleosl('demo1',xver,th,p_m1,[],filesuff);
            thhat_m1=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    try
        thhat_m1_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_m1_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_m1,[],filesuff,th,[1 0 1]);
        thhat_m1_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_m1_f,1)<xver 
            mleosl('demo1',xver,th,p_m1,[],filesuff,th,[1 0 1]);
            thhat_m1_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The mask exterior case
    p_m2=p;p_m2.taper=~I;
    filesuff=append('mask2_',datum);
    try
        thhat_m2=load(append('mleosl_thhat_',filesuff));
        if size(thhat_m2,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_m2,[],filesuff);
        thhat_m2=load(append('mleosl_thhat_',filesuff));
        while size(thhat_m2,1)<xver 
            mleosl('demo1',xver,th,p_m2,[],filesuff);
            thhat_m2=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    try
        thhat_m2_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_m2_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_m2,[],filesuff,th,[1 0 1]);
        thhat_m2_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_m2_f,1)<xver 
            mleosl('demo1',xver,th,p_m2,[],filesuff,th,[1 0 1]);
            thhat_m2_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The speckle case
    p_s1=p;p_s1.mask='random';
    [~,I]=muckit(ones(p_s1.NyNx),p_s1);
    p_s1.taper=I;
    filesuff=append('muck1_',datum);
    try
        thhat_s1=load(append('mleosl_thhat_',filesuff));
        if size(thhat_s1,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_s1,[],filesuff);
        thhat_s1=load(append('mleosl_thhat_',filesuff));
        while size(thhat_s1,1)<xver 
            mleosl('demo1',xver,th,p_s1,[],filesuff);
            thhat_s1=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    try
        thhat_s1_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_s1_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_s1,[],filesuff,th,[1 0 1]);
        thhat_s1_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_s1_f,1)<xver 
            mleosl('demo1',xver,th,p_s1,[],filesuff,th,[1 0 1]);
            thhat_s1_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The anti-speckle case
    p_s2=p;p_s2.mask='random';p_s2.taper=~I;
    filesuff=append('muck2_',datum);
    try
        thhat_s2=load(append('mleosl_thhat_',filesuff));
        if size(thhat_s2,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_s2,[],filesuff);
        thhat_s2=load(append('mleosl_thhat_',filesuff));
        while size(thhat_s2,1)<xver 
            mleosl('demo1',xver,th,p_s2,[],filesuff);
            thhat_s2=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    try
        thhat_s2_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_s2_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_s2,[],filesuff,th,[1 0 1]);
        thhat_s2_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_s2_f,1)<xver 
            mleosl('demo1',xver,th,p_s2,[],filesuff,th,[1 0 1]);
            thhat_s2_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The setup for creating the plots of the simulated fields is very similar
    % to est==0, but we will need to load our previously stored  XVER estimates
    % from their MLEOSL hash
    
    % maxsd needs to be calculated for every set of thhats:  
    maxsd=max(max([std(thhat_r1);std(thhat_r1_f);std(thhat_r2);std(thhat_r2_f);
               std(thhat_m1);std(thhat_m1_f);std(thhat_m2);std(thhat_m2_f);
               std(thhat_s1);std(thhat_s1_f);std(thhat_s2);std(thhat_s2_f)]));

    % Figure formatting
    % Font size
    afs=10;
    % Marker size
    ms=100;
    
    figure(1)
    clf
    [ah,ha,H]=krijetem(subnum(4,3));
    Hx=simulosl(th,p);
    for hnd=1:3
        if hnd==1
            % The rectangular case; for the full field, plot vw; for the
            % external region of the rectangle, plot w
            [v,I,w,vw,scl,cr]=maskit(Hx,p_r2,scl,Hx);
            top=vw;bottom=w;mt=0;
            thhat1=thhat_r1;thhat1f=thhat_r1_f;thhat2=thhat_r2;thhat2f=thhat_r2_f;
        elseif hnd==2
            % The geographic boundary case; 
            [v,I,w,vw,scl,cr]=maskit(Hx,p_m1,scl,Hx);
            top=w;bottom=Hx;bottom(~~I(:))=NaN;mt=0;
            thhat1=thhat_m1;thhat1f=thhat_m1_f;thhat2=thhat_m2;thhat2f=thhat_m2_f;
        elseif hnd==3
            % The speckle case
            [v,I,w,scl,cr]=muckit(Hx,p_s1,scl);
            top=v;bottom=w;mt=1;
            thhat1=thhat_s1;thhat1f=thhat_s1_f;thhat2=thhat_s2;thhat2f=thhat_s2_f;
        else                                                                    
            error('there is an indexing error; debug creation of subfigures')   
        end                                                                     
        axes(ha((hnd-1)*4+1))                                                   
        t=plotit(top,p,cr,th,mt,afs);
        %if ~isempty(t)
        %    t.Position(1)=-3*t.Position(1);
            %t.Position(1)=-1.75*t.Position(1);
        %end

        axes(ha((hnd-1)*4+2))
        % Scatter marker for each th0 value, scale will be standardized by +/-3
        % std dev
        symbs={'^','s','o'};
        yvls=[5 3 1;nan(2,3)];yfvls=[5 nan 1;nan(2,3)];nxvls=nan(2,3);
        hold on
        % th0 values:
        th0p=scatter([0 0 0;nxvls],yvls,'k','filled');
        for ind=1:3;th0p(ind).Marker=symbs{ind};end
        % standardize estimates:
        disp(sprintf('only considering the first %i estimates, additional estimates will be ignored',xver))
        t1m=mean(thhat1(1:xver,:)-th);t1fm=mean(thhat1f(1:xver,:)-th);
        t1sd=std(thhat1(1:xver,:));t1fsd=std(thhat1f(1:xver,:));
       % maxsd needs to be calculated for every set of thhats: maxsd=max([t1sd,t1fsd]);
        t1sd=t1sd./maxsd;t1fsd=t1fsd./maxsd;
        maxdst=max(abs([t1m+t1sd,t1m-t1sd,t1fm+t1fsd,t1fm-t1fsd]));
        thhat1fsp=plot([t1fm-t1fsd;t1fm+t1fsd],[yfvls(1,:);yfvls(1,:)],'b-');
        thhat1sp=plot([t1m-t1sd;t1m+t1sd],[yvls(1,:);yvls(1,:)],'k-');
        thhat1fspb=plot([t1fm-t1fsd;t1fm+t1fsd],[yfvls(1,:);yfvls(1,:)],'b|');
        thhat1spb=plot([t1m-t1sd;t1m+t1sd],[yvls(1,:);yvls(1,:)],'k|');
        % Let's look at creating kde plots to show the shape of the
        % distribution, assuming a Gaussian kernel:
        plotkde=1;
        if plotkde==1
             [bandwidth,density,xmesh]=kde(thhat1(1:xver,1)-th(1));
             plot(xmesh,density+5,'k')
             [bandwidth,density,xmesh]=kde(thhat1(1:xver,2)-th(2));
             plot(xmesh,density+3,'k')
             [bandwidth,density,xmesh]=kde(thhat1(1:xver,3)-th(3));
             plot(xmesh,density+1,'k')
             [bandwidth,density,xmesh]=kde(thhat1f(1:xver,1)-th(1));
             plot(xmesh,density+5,'b')
             [bandwidth,density,xmesh]=kde(thhat1f(1:xver,3)-th(3));
             plot(xmesh,density+1,'b')
        elseif plotkde==-1
              try                                                                           
                  [a,bdens,c]=kde(thhat1(:,ind),2^8);                                       
              catch                                                                         
                  a=Inf;                                                                    
              end                                                                           
              if isinf(a) || any(bdens<-1e-10) || size(thhat1,1)<50                         
                [bdens,c]=hist(thhat1(:,ind),max(size(thhat1,1)/3,3));                      
                bdens=bdens/indeks(diff(c),1)/size(thhat1(:,ind),1);                        
              end                                                                           
              % This number is close to one... it's a proper density!                       
              if xver==2                                                                    
                disp(sprintf('%s pdf normalization check by summation %g',...               
                             upper(mfilename),sum(bdens)*indeks(diff(c),1)))                
                disp(' ')                                                                   
              end                                                                           
                                                                                            
              % Note that there are changes to plotting functions in 2023b
              % Plot using a scale factor to remove the units from the y axis        
              thhis(ind)=bar(c,t1sd*bdens,1);                                               
              %fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C );
              %      set(ah(ind),'ylim',yls)                                                      
              % The markings                                                                
              %stats=mobs+nstats*t1sd;                                                       
              % What is the percentage captured within the range?                           
              %nrx=20; nry=15;                                                               
              %try                                                                           
              %   set(ah(ind),'XLim',stats([1 end]),'XTick',stats,'XTickLabel',nstats)         
              %end 
        end

        %thhat1p=scatter([t1m;nxvls],yvls,'k','filled');
        %thhat1fp=scatter([t1fm;nxvls],yvls,'r','filled');
        %for ind=1:3;thhat1p(ind).Marker=symbs{ind};thhat1fp(ind).Marker=symbs{ind};end
        
        xlim([-maxdst,maxdst]*1.1);
        ylim([1,6]+1*[-1,1]);
        h=gca;h.XAxis.Visible='off';h.YAxis.Visible='off';
        %text(-maxdst/2,6,sprintf('$\\sigma^2_0=%0.1f$',th(1)),'Interpreter', 'latex');
        %text(-maxdst/2,4,sprintf('$\\nu_0=%0.1f$',th(2)),'Interpreter', 'latex');
        %text(-maxdst/2,2,sprintf('$\\rho_0=%0.1f$',th(3)),'Interpreter', 'latex');
        hold off

        axes(ha((hnd-1)*4+3))                                                   
        plotit(bottom,p,cr,th,0,afs);

        axes(ha((hnd-1)*4+4))
        hold on
        % th0 values:
        th0p=scatter([0 0 0;nxvls],yvls,'k','filled');
        for ind=1:3;th0p(ind).Marker=symbs{ind};end
        t2m=mean(thhat2(1:xver,:)-th);t2fm=mean(thhat2f(1:xver,:)-th);
        t2sd=std(thhat2(1:xver,:));t2fsd=std(thhat2f(1:xver,:));
        % maxsd needs to be calculated for every set of thhats:  maxsd=max([t2sd,t2fsd]);
        t2sd=t2sd./maxsd;t2fsd=t2fsd./maxsd;
        maxdst=max(abs([t2m+t2sd,t2m-t2sd,t2fm+t2fsd,t2fm-t2fsd]));
        thhat2fsp=plot([t2fm-t2fsd;t2fm+t2fsd],[yfvls(1,:);yfvls(1,:)],'b-');
        thhat2sp=plot([t2m-t2sd;t2m+t2sd],[yvls(1,:);yvls(1,:)],'k-');
        thhat2fspb=plot([t2fm-t2fsd;t2fm+t2fsd],[yfvls(1,:);yfvls(1,:)],'b|');
        thhat2spb=plot([t2m-t2sd;t2m+t2sd],[yvls(1,:);yvls(1,:)],'k|');
        %thhat1p=scatter([t1m;nxvls],yvls,'k','filled');
        %thhat1fp=scatter([t1fm;nxvls],yvls,'r','filled');
        %for ind=1:3;thhat1p(ind).Marker=symbs{ind};thhat1fp(ind).Marker=symbs{ind};end

        if plotkde==1
            [bandwidth,density,xmesh]=kde(thhat2(1:xver,1)-th(1));
            plot(xmesh,density+5,'k')
            [bandwidth,density,xmesh]=kde(thhat2(1:xver,2)-th(2));
            plot(xmesh,density+3,'k')
            [bandwidth,density,xmesh]=kde(thhat2(1:xver,3)-th(3));
            plot(xmesh,density+1,'k')
            [bandwidth,density,xmesh]=kde(thhat2f(1:xver,1)-th(1));
            plot(xmesh,density+5,'b')
            [bandwidth,density,xmesh]=kde(thhat2f(1:xver,3)-th(3));
            plot(xmesh,density+1,'b')
        elseif plotkde==-1
              try                                                                           
                  [a,bdens,c]=kde(thhats(:,ind),2^8);                                       
              catch                                                                         
                  a=Inf;                                                                    
              end                                                                           
              if isinf(a) || any(bdens<-1e-10) || size(thhats,1)<50                         
                % If it's funky use the good old histogram                                  
                [bdens,c]=hist(thhats(:,ind),max(size(thhats,1)/3,3));                      
                bdens=bdens/indeks(diff(c),1)/size(thhats(:,ind),1);                        
              end                                                                           
              % This number is close to one... it's a proper density!                       
              if xver==2                                                                    
                disp(sprintf('%s pdf normalization check by summation %g',...               
                             upper(mfilename),sum(bdens)*indeks(diff(c),1)))                
                disp(' ')                                                                   
              end                                                                           
                                                                                            
              % Note that there are changes that wreck this in 2023b (?)                    
                                                                                            
              % Now plot it using a scale factor to remove the units from the y axis        
              thhis(ind)=bar(c,sobs*bdens,1);                                               
              % set(ah(ind),'ylim',yls)                                                       
              % The markings                                                                
              stats=mobs+nstats*sobs;                                                       
              % What is the percentage captured within the range?                           
              nrx=20; nry=15;                                                               
              % try                                                                           
              %     set(ah(ind),'XLim',stats([1 end]),'XTick',stats,'XTickLabel',nstats)         
              % end 
        end

        xlim([-maxdst,maxdst]*1.1);
        ylim([1,6]+1*[-1,1]);
        h=gca;h.XAxis.Visible='off';h.YAxis.Visible='off';
        %text(-maxdst/2,6,sprintf('$\\sigma^2_0=%0.1f$',th(1)),'Interpreter', 'latex');
        %text(-maxdst/2,4,sprintf('$\\nu_0=%0.1f$',th(2)),'Interpreter', 'latex');
        %text(-maxdst/2,2,sprintf('$\\rho_0=%0.1f$',th(3)),'Interpreter', 'latex');
        hold off
    end  
    figdisp([],sprintf('samplingfig2_%3.3i_%3.3i_%3.3i_fld',p.NyNx,xver),[],1);
elseif est==2
    % Number of iterations
    defval('xver',100);
    datum='18-Aug-2025';%date;

    % This is a nice grid size for visuals 
    p=[];p.NyNx=[188 233];p.dydx=[5e3 5e3];p.blurs=Inf;
    % This NyNx matches structured mask gebcostruct
    p=[];p.NyNx=[117 124];p.dydx=[5e3 5e3];p.blurs=Inf;
    % This NyNx matches structured mask gebcostruct0 following decimation x5
    p=[];p.NyNx=[261 255];p.dydx=[5e3 5e3];p.blurs=Inf;
    th=[1e3 1 2e4];

    % 0 correct, 1 incorrect fix
    defval('fixflag',0)
    thfi=[th(1) 1/2 th(3)];
    % default correct incorrect
    thfs={[] th thfi};
    % 3p 2p 2p
    ifinvs={[1 1 1] [1 0 1] [1 0 1]};

    % We will shuffle through three tapers and their complements, and we will
    % have three variable name styles for each (3p,2pc,2pi)
    fsuf={sprintf('bound_3p_%s',datum),sprintf('struct_3p_%s',datum),sprintf('rand_3p_%s',datum),...
          sprintf('abound_3p_%s',datum),sprintf('astruct_3p_%s',datum),sprintf('arand_3p_%s',datum),...
          sprintf('bound_2c_%s',datum),sprintf('struct_2c_%s',datum),sprintf('rand_2c_%s',datum),...
          sprintf('abound_2c_%s',datum),sprintf('astruct_2c_%s',datum),sprintf('arand_2c_%s',datum),...
          sprintf('bound_2i_%s',datum),sprintf('struct_2i_%s',datum),sprintf('rand_2i_%s',datum),...
          sprintf('abound_2i_%s',datum),sprintf('astruct_2i_%s',datum),sprintf('arand_2i_%s',datum)};
    ps={p p p p p p};
    ps{1}.mask='newjersey';
    [~,I]=maskit(rand(p.NyNx),ps{1},0.662);
    ps{1}.taper=double(I);
    ps{2}.mask='antinewjersey';
    ps{2}.taper=double(~I);

    ps{3}.mask='structured';
    % 0.336 observed, NyNx=[261 255]
    geb=load('MYMASKS/gebcostruct0.mat');
    decm=5;
    I=geb.gebcostruct0(1:decm:end,1:decm:end);
    ps{3}.taper=double(I);
    ps{4}.mask='antistructured';
    ps{4}.taper=double(~logical(I)); 

    ps{5}.mask='random';
    [~,I]=muckit(rand(p.NyNx),ps{5},0.336);
    ps{5}.taper=double(I);
    ps{6}.mask='antirandom';
    ps{6}.taper=double(~I);

    % We can also perform some calculations for the full rectangle case, but we
    % will not show these results in the figure
    pr=p;pr.taper=1;
    fsufr={sprintf('full_3p_%s',datum),sprintf('full_2c_%s',datum),sprintf('full_2i_%s',datum)};
    thhatr=cell(3,1);
    for ind=1:3
      try
        thhatr{ind}=load(sprintf('mleosl_thhat_%s',fsufr{ind}));
        if size(thhatr{ind},1)<xver;error('Calculating additional realizations');end
      catch
        cxver=xver-size(thhatr{ind},1);
        mleosl('demo1',cxver,th,pr,[],fsufr{ind},...
               thfs{ind},ifinvs{ind});
        thhatr{ind}=load(sprintf('mleosl_thhat_%s',fsufr{ind}));
        while size(thhatr{ind},1)<xver 
          cxver=xver-size(thhatr{ind},1);
          mleosl('demo1',cxver,th,pr,[],fsufr{ind},...
                 thfs{ind},ifinvs{ind});
          thhatr{ind}=load(sprintf('mleosl_thhat_%s',fsufr{ind}));
        end
      end
      stdobsr{ind}=std(thhatr{ind})';
      mobsr{ind}=mean(thhatr{ind})';
      % Calculate covariance for the truth, or assumed truth in the incorrect
      % special case case, with the current mask
      cm=1;
      if isempty(thfs{ind})
        cth=th;
      else
        cth=thfs{ind};
      end
      pr.blurs=-1;
      covt=covthosl(cth,pr,cm,ifinvs{ind});
      covth{ind}=matslice(covt,ifinvs{ind},-1);
      covsd{ind}=sqrt(diag(covth{ind}));
      pr.blurs=Inf;
    end

    thhats=cell(18,1);
    % 6 sampling grids, 3 covariance parameter settings; work through all
    % sampling grids for 3p, then all grids for 2c, and finally all grids for 2i
    for ind=1:6*3
      try
          thhats{ind}=load(sprintf('mleosl_thhat_%s',fsuf{ind}));
          if size(thhats{ind},1)<xver;error('Calculating additional realizations');end
      catch
          cxver=xver-size(thhats{ind},1);
          mleosl('demo1',cxver,th,ps{mod(ind-1,6)+1},[],fsuf{ind},...
                 thfs{mod(floor((ind-1)/6),3)+1},ifinvs{mod(floor((ind-1)/6),3)+1});
          thhats{ind}=load(append('mleosl_thhat_',fsuf{ind}));
          while size(thhats{ind},1)<xver 
            cxver=xver-size(thhats{ind},1);
            mleosl('demo1',cxver,th,ps{mod(ind-1,6)+1},[],fsuf{ind},...
              thfs{mod(floor((ind-1)/6),3)+1},ifinvs{mod(floor((ind-1)/6),3)+1});
            thhats{ind}=load(append('mleosl_thhat_',fsuf{ind}));
          end
      end
      stdobs{ind}=std(thhats{ind})';
      mobs{ind}=mean(thhats{ind})';
      % Calculate covariance for the truth, or assumed truth in the incorrect
      % special case case, with the current mask
      cm=1;
      if isempty(thfs{mod(floor((ind-1)/6),3)+1})
        cth=th;
      else
        cth=thfs{mod(floor((ind-1)/6),3)+1};
      end
      ps{mod(ind-1,6)+1}.blurs=-1;
      covt=covthosl(cth,ps{mod(ind-1,6)+1},cm,ifinvs{mod(floor((ind-1)/6),3)+1});
      covth{ind}=matslice(covt,ifinvs{mod(floor((ind-1)/6),3)+1},-1);
      covsd{ind}=sqrt(diag(covth{ind}));
      ps{mod(ind-1,6)+1}.blurs=Inf;
    end
    thfs{1}=th;

    % for the table:
    % tnd=2;
    % mobs{tnd}(1), covsd{tnd}(1)/th(1)*100, stdobs{tnd}(1)/th(1)*100
    % mobs{tnd}(2), covsd{tnd}(2)/th(2)*100, stdobs{tnd}(2)/th(2)*100
    % pi*mobs{tnd}(3)/L, covsd{tnd}(3)/th(3)*100, stdobs{tnd}(3)/th(3)*100
    % covth{tnd}./[diag(sqrt(covth{tnd}))*diag(sqrt(covth{tnd}))']*100
    % covob{tnd}=cov(thhats{tnd});
    % covob{tnd}./[diag(sqrt(covob{tnd}))*diag(sqrt(covob{tnd}))']*100

    % largest parameter spread
    maxsd=max(max([covsd{:}]'),max([stdobs{:}]'));
    datspc=0.1;
    stdOut=4;
    % Labels for the nu axes will be based on the special values
    nutickvals=[1/3 1/2 1 3/2 5/2];
    %nuticklabs={'$1/3$','$1/2$','$1$','$\nu_0$','$5/2$'};
    nuticklabs={'$\frac{1}{3}$','$\frac{1}{2}$','$1$','$\frac{3}{2}$','$\frac{5}{2}$'};
    triy=-0.06;
    L=floor(sqrt(sum((p.NyNx.*p.dydx).^2)));
    % color options for kde plots
    clrs=["#808080","#7570B3","#EDB120","#0072BD","#1B9E77"];
    % face and edge alpha
    fa=0.25;
    % font size
    afs=10;
    % marker size
    ms=100;
    % power of 2 for kde
    kdep=9;
    figure(2)
    clf
    numrows=8;
    labs={'$\hat{\sigma}^2$','$\hat{\nu}$','$\hat{\rho}$'};
    [ah,ha,H]=krijetem(subnum(numrows,3));
    Hx=simulosl(th,p);
    % 6 sampling grids 
    for jnd=1:6
      cf=Hx.*ps{jnd}.taper(:);
      thhat3p=thhats{jnd};
      covsd3p=covsd{jnd};
      mobs3p=mobs{jnd};
      obssd3p=stdobs{jnd};
      thhat2c=thhats{jnd+6};
      covsd2c=covsd{jnd+6};
      mobs2c=mobs{jnd+6};
      obssd2c=stdobs{jnd+6};
      thhat2i=thhats{jnd+12};
      covsd2i=covsd{jnd+12};
      mobs2i=mobs{jnd+12};
      obssd2i=stdobs{jnd+12};
      %
      axes(ha((jnd-1)*numrows/2+1)) 
      hold on
      % cr comes out of maskit, can we get around this?
      cr=[];
      if ismember(jnd,[3 4]);mt=1;else;mt=0;end
      [h,t]=plotit(cf,ps{jnd},cr,[],mt,afs);
      longticks
      box on
      if jnd>1;set(gca,'YTickLabel',[]);end
      hold off
      % 3 parameters
      for knd=1:3
        axes(ha((jnd-1)*numrows/2+1+knd))
        hold on
        datspc=th(knd);%/2;
        % maybe we want everything on the same axis
        kdelim3{jnd,knd}=[min(thhat3p(:,knd)) max(thhat3p(:,knd))]+[-1 1].*datspc;
        kdelim3{jnd,knd}=th(knd)+[-maxsd(:,knd)*stdOut maxsd(:,knd)*stdOut];
        % kdelimc{jnd,knd}=[min(thhat2c(:,knd)) max(thhat2c(:,knd))]+[-1 1].*datspc;
        % kdelimi{jnd,knd}=[min(thhat2i(:,knd)) max(thhat2i(:,knd))]+[-1 1].*datspc;
        [a113,bdens113,c113]=kde(thhat3p(:,knd),2^kdep,kdelim3{jnd,knd}(1),kdelim3{jnd,knd}(2));
        % Scale the kde to remove units from the vertical axis
        cnd=3*(jnd-1)+1;
        %t113kp=bar(c113,bdens113*stdobs{cnd}(knd),1,'FaceColor',clrs{1},...
        t113kp=bar(c113,bdens113*obssd3p(knd),1,'FaceColor',clrs{1},...
                     'EdgeColor',clrs{1},'FaceAlpha',fa,'EdgeAlpha',fa);
        % Observed and predicted covariance curves for estimated x-axis; follows
        % from MLEPLOS
        % xnorm3p=linspace(-stdOut,stdOut,100).*stdobs{cnd}(knd)+mobs{cnd}(knd);
        xnorm3p=linspace(kdelim3{jnd,knd}(1),kdelim3{jnd,knd}(2),300);
        xnorm3p=th(knd)+linspace(-maxsd(:,knd)*stdOut*3,maxsd(:,knd)*stdOut*3,500);
        % pobs3p{jnd,knd}=plot(xnorm3p,stdobs{cnd}(knd).*...
        %                 normpdf(xnorm3p,mobs{cnd}(knd),stdobs{cnd}(knd)));
        pobs3p{jnd,knd}=plot(xnorm3p,obssd3p(knd).*...
                        normpdf(xnorm3p,mobs3p(knd),obssd3p(knd)));
        % pclc3p{jnd,knd}=plot(xnorm3p,covsd{cnd}(knd).*...
        %                 normpdf(xnorm3p,mobs{cnd}(knd),covsd{cnd}(knd)));
        pclc3p{jnd,knd}=plot(xnorm3p,covsd3p(knd).*...
                        normpdf(xnorm3p,mobs3p(knd),covsd3p(knd)));
        set(pclc3p{jnd,knd},'LineWidth',0.5,'Color','k')
        set(pobs3p{jnd,knd},'LineWidth',1.0,'Color',grey(3.5))
        if knd~=2
          [a11c,bdens11c,c11c]=kde(thhat2c(:,knd),2^kdep,kdelim3{jnd,knd}(1),kdelim3{jnd,knd}(2));
          cnd=3*(jnd-1)+2;
          %t11ckp=bar(c11c,-1*bdens11c*stdobs{cnd}(knd),1,'FaceColor',clrs{2},...
          t11ckp=bar(c11c,-1*bdens11c*obssd2c(knd),1,'FaceColor',clrs{2},...
                     'EdgeColor',clrs{2},'FaceAlpha',fa,'EdgeAlpha',fa);
          % xnorm2c=linspace(-stdOut,stdOut,100).*stdobs{cnd}(knd)+mobs{cnd}(knd);
          % pobs2c{jnd,knd}=plot(xnorm3p,-stdobs{cnd}(knd).*...
          %                 normpdf(xnorm3p,mobs{cnd}(knd),stdobs{cnd}(knd)));
          % pclc2c{jnd,knd}=plot(xnorm3p,-covsd{cnd}(knd).*...
          %                 normpdf(xnorm3p,mobs{cnd}(knd),covsd{cnd}(knd)));
          pobs2c{jnd,knd}=plot(xnorm3p,-obssd2c(knd).*...
                          normpdf(xnorm3p,mobs2c(knd),obssd2c(knd)));
          pclc2c{jnd,knd}=plot(xnorm3p,-covsd2c(knd).*...
                          normpdf(xnorm3p,mobs2c(knd),covsd2c(knd)));
          set(pclc2c{jnd,knd},'LineWidth',0.5,'Color','k')
          set(pobs2c{jnd,knd},'LineWidth',1.0,'Color',clrs{2})
          %
          [a11i,bdens11i,c11i]=kde(thhat2i(:,knd),2^kdep,kdelim3{jnd,knd}(1),kdelim3{jnd,knd}(2));
          cnd=3*(jnd-1)+3;
          %t11ikp=bar(c11i,-1*bdens11i*stdobs{cnd}(knd),1,'FaceColor',clrs{3},...
          t11ikp=bar(c11i,-1*bdens11i*obssd2i(knd),1,'FaceColor',clrs{3},...
                     'EdgeColor',clrs{3},'FaceAlpha',fa,'EdgeAlpha',fa);
          % xnorm2i=linspace(-stdOut,stdOut,100).*stdobs{cnd}(knd)+mobs{cnd}(knd);
          % pobs2i{jnd,knd}=plot(xnorm3p,-stdobs{cnd}(knd).*...
          %                 normpdf(xnorm3p,mobs{cnd}(knd),stdobs{cnd}(knd)));
          % pclc2i{jnd,knd}=plot(xnorm3p,-covsd{cnd}(knd).*...
          %                 normpdf(xnorm3p,mobs{cnd}(knd),covsd{cnd}(knd)));
          pobs2i{jnd,knd}=plot(xnorm3p,-obssd2i(knd).*...
                          normpdf(xnorm3p,mobs2i(knd),obssd2i(knd)));
          %Does it make sense to look at the predicted covariance calculated for
          %an incorrect model...? 
          % pclc2i{jnd,knd}=plot(xnorm3p,-covsd2i(knd).*...
          %                 normpdf(xnorm3p,mobs2i(knd),covsd2i(knd)));
          set(pobs2i{jnd,knd},'LineWidth',1.0,'Color',clrs{3})
         %set(pclc2i{jnd,knd},'LineWidth',0.5,'Color','k')
        end
  
      if knd==2
          th02fc{jnd}=plot(thfs{2}(2),triy,'Marker','v',...
                      'MarkerFaceColor',clrs{2},'MarkerEdgeColor','k',...
                      'MarkerSize',8);
          th02fi{jnd}=plot(thfs{3}(2),triy,'Marker','v',...
                      'MarkerFaceColor',clrs{3},'MarkerEdgeColor','k',...
                      'MarkerSize',8);
          maxy(jnd,knd)=max(abs([bdens113*stdobs{cnd}(knd)]));
          maxx(jnd,knd)=max(abs([c113]));
        else
          maxx(jnd,knd)=nanmax(abs([c113 c11c c11i]));
          %maxy(jnd,knd)=max(abs([bdens113*stdobs{3*(jnd-1)+1}(knd);...
          %                       bdens11c*stdobs{3*(jnd-1)+2}(knd);...
          %                       bdens11i*stdobs{3*(jnd-1)+3}(knd)]));
          maxy(jnd,knd)=nanmax(abs([bdens113*obssd3p(knd);...
                                    bdens11c*obssd2c(knd);...
                                    bdens11i*obssd2i(knd)]));
        end
        th013{jnd,knd}=scatter(th(knd),0,ms,'k|');
        %plot([min(abs([c113 c11c c11i])),maxx(jnd,knd)],[0,0],'k')
        plot([mean([c113 c11c c11i])+stdOut.*maxsd(:,knd).*[-1 1]],[0 0],'k')
        if jnd<3
          ylabel(labs{knd},'Interpreter', 'latex')
          set(get(gca,'ylabel'),'rotation',0)
        end
        longticks
        set(gca,'YTickLabel',[])
        box on
        hold off
  
        if knd==2
          set(gca,'XTick',nutickvals,'XTickLabel',nuticklabs,...
              'TickLabelInterpreter','latex','FontSize',afs,...
              'XTickLabelRotationMode','manual');
        end
        box on
        hold off
      end
    end
    % Now assign a common axis range for all of the kde subplots
    for knd=1:3
        if knd==2
          ylspc=0.5;
        else
          ylspc=1.1;
        end
        %set(H(hnd+1,:),'XLim',th(hnd)+(max(maxx(:,hnd))-th(hnd)).*[-1*xlspc xlspc],...
        set(H(knd+1,:),'XLim',th(knd)+maxsd(:,knd).*[-1*stdOut stdOut],...
          'YLim',max(maxy(:,knd)).*[-1*ylspc ylspc]);
        set(H(knd+5,:),'XLim',th(knd)+maxsd(:,knd).*[-1*stdOut stdOut],...
          'YLim',max(maxy(:,knd)).*[-1*ylspc ylspc]);
    end

    % Treat labels for sigma and rho
    for hnd=1:3
        for jnd=1:2
        % select the sigma plots
        axes(H(2+4*(jnd-1),hnd));
        xtickvals=gca().XAxis.TickValues;
        if any(xtickvals==th(1));
            xticklabs=cellstr(num2str(xtickvals(:),'%0.1g'));
            xticklabs{xtickvals==th(1)}=sprintf('$\\sigma^2_0$=%0.1g',th(1));
            xticklabs={xticklabs{1}; xticklabs{xtickvals==th(1)};...
              xticklabs{end}};
            xtickvals=[xtickvals(1) th(1) xtickvals(end)];
        else
            [xtickvals,srtsig]=sort([xtickvals th(1)]);
            xticklabs=cellstr(num2str(xtickvals(:),'%0.1g'));
            xticklabs{srtsig==max(srtsig)}=sprintf('$\\sigma^2_0$=%0.1g',th(1));
            xticklabs={xticklabs{1}; xticklabs{srtsig==max(srtsig)};...
              xticklabs{end}};
            xtickvals=[xtickvals(1) th(1) xtickvals(end)];
        end
        set(gca,'XTick',xtickvals,'XTickLabel',xticklabs,...
            'TickLabelInterpreter','latex','FontSize',afs,...
            'XTickLabelRotationMode','manual');

        % select the rho plots
        axes(H(4+4*(jnd-1),hnd));
        xtickvals=gca().XAxis.TickValues;
        if any(xtickvals==th(3));
            xticklabs=cellstr(num2str(xtickvals(:)./L,'%0.1g'));
            xticklabs{xtickvals==th(3)}=sprintf('$\\rho_0/L$=%0.1g',th(3)/L);
            xticklabs={xticklabs{1}; xticklabs{xtickvals==th(3)};...
                xticklabs{end}};
            xtickvals=[xtickvals(1) th(3) xtickvals(end)];
        else
          [xtickvals,srtrh]=sort([xtickvals th(3)]);
          xticklabs=cellstr(num2str(xtickvals(:)./L,'%0.1g'));
          xticklabs{srtrh==max(srtrh)}=sprintf('$\\rho_0/L$=%0.1g',th(3)/L);
          xticklabs={xticklabs{1}; xticklabs{srtrh==max(srtrh)};...
              xticklabs{end}};
          xtickvals=[xtickvals(1) th(3) xtickvals(end)];
        end
        set(gca,'XTick',xtickvals,'XTickLabel',xticklabs,...
          'TickLabelInterpreter','latex','FontSize',afs,...
          'XTickLabelRotationMode','manual');
    end

        % Special treatment for nu x-axis
        xlspc=1.1;
        set(H(3,:),'XLim',th(2)+(max(abs(th(2)-nutickvals))).*[-1*xlspc xlspc]);
        set(H(7,:),'XLim',th(2)+(max(abs(th(2)-nutickvals))).*[-1*xlspc xlspc]);
    end  

    axes(H(1,1));rf1pos=gca().Position;axes(H(1,2));rf2pos=gca().Position;
    axes(H(1,3));rf3pos=gca().Position; 
    axes(H(5,1));rf4pos=gca().Position;axes(H(5,2));rf5pos=gca().Position;
    axes(H(5,3));rf6pos=gca().Position;
    fshr=0.5;
    rf1npos=shrink(H(1,1),fshr,fshr);
    rf2npos=shrink(H(1,2),fshr,fshr);
    rf3npos=shrink(H(1,3),fshr,fshr);
    % For some reason SHRINK isn't always effective, so the below lines are
    % needed for the latter fields to have an equivalent size to the first
    for hnd=2:3
         ha((hnd-1)*numrows+1).Position(3)=ha(1).Position(3);
         ha((hnd-1)*numrows+1).Position(4)=ha(1).Position(4);
         ha((hnd-1)*numrows+1).Position(2)=ha(1).Position(2);
    end

    H(5,1).Position(3)=ha(1).Position(3);
    H(5,1).Position(4)=ha(1).Position(4);
    H(5,1).Position(1)=ha(1).Position(1);
    for hnd=2:3
         H(5,hnd).Position(3)=H(5,1).Position(3);
         H(5,hnd).Position(4)=H(5,1).Position(4);
         H(5,hnd).Position(2)=H(5,1).Position(2);
    end

    % move the rows of fields down to their original alignment
    H(1,1).Position(2)=rf1pos(2)-rf1npos(4)+rf1pos(4);
    H(1,2).Position(2)=H(1,1).Position(2);
    H(1,3).Position(2)=H(1,1).Position(2);
    H(5,1).Position(2)=rf4pos(2)-rf1npos(4)+rf4pos(4);
    H(5,2).Position(2)=H(5,1).Position(2);
    H(5,3).Position(2)=H(5,1).Position(2);
    % move the bottom right fields to be centered with the top fields
    H(5,2).Position(1)=H(1,2).Position(1);
    H(5,3).Position(1)=H(1,3).Position(1);

    thshr=1.5;
    shrink([H(2:4,:)],1,thshr)
    serre([H(2:4,1)],1/4,'down');serre([H(2:4,2)],1/4,'down');serre([H(2:4,3)],1/4,'down')
    movev(H(1,:),0.05);
    shrink([H(6:8,:)],1,thshr)
    serre([H(6:8,1)],1/4,'down');serre([H(6:8,2)],1/4,'down');serre([H(6:8,3)],1/4,'down')
    movev(H(6:8,:),-0.05);
    movev(H,-0.04);

    % Minimize space between the subplots
    serre([H(2:4,:)],1/2,'across')
    serre([H(6:8,:)],1/2,'across')
    serre([H(1,:)'],-15/64,'across')
    serre([H(5,:)'],-15/64,'across')
    movev(H(1:4,:),0.01)
    H(1,1).XTickLabelRotation=0;H(1,2).XTickLabelRotation=0;H(1,3).XTickLabelRotation=0;
    H(5,1).XTickLabelRotation=0;H(5,2).XTickLabelRotation=0;H(5,3).XTickLabelRotation=0;
    H(5,1).YTickLabel=H(1,1).YTickLabel;
    axes(H(1,1));box on;axes(H(1,2));box on;axes(H(1,3));box on;
    axes(H(5,1));box on;axes(H(5,2));box on;axes(H(5,3));box on;

    % Make some adjustments for revision
    for ind=[1 2 3 13 14 15]
      ah(ind).XTickLabels=[1 64 127.5];
      ah(ind).XLabel.String='easting [km]';
      ah(ind).FontSize=8;
      if ind==1|ind==13
        ah(ind).YTickLabels=[1 66 130.5];
        ah(ind).YLabel.String='northing [km]';
      end
    end
    axes(ah(3))
    [cb,xcb]=addcb('vert',halverange(Hx,80,NaN),halverange(Hx,80,NaN),'parula',[],0);
    xcb.String='elevation [m]';
    cb.TickDir='out';     
    cb.YAxisLocation='right';
    cb.YTick=linspace(cb.YTick(1),cb.YTick(end),7);
    hvr=halverange(Hx,80,NaN);
    cb.YTickLabel=round(linspace(hvr(1),hvr(2),7),0);
    cb.Position(1)=0.89;
    cb.Position(3)=0.01;

    [bh,txh]=label([ah(1:3) ah(13:15)],'ul',8,0);
    set(gcf,'Renderer','painters')

    fname=sprintf('samplingfig2_%3.3i_%3.3i_%3.3i_fld_v3',p.NyNx,xver);
    figdisp([],fname,[],1);
elseif est==3
    % Figure comparing full rectangle to inner an inner-rectangle padded by
    % zeros and a subfield of only the inner-rectangle

    % Once we find a grid-size we like through experimentation, let's set this
    % size for future program calls
    p.NyNx = [188, 233]; 
    th = [1 1.5 5];

    figure(3)
    clf
    numrows=4;
    [ah,ha,H]=krijetem(subnum(numrows,3));
    %Hx=simulosl(th,p);
    Hx=simulosl(th,p);
    p_r1=p; % Save the original taper for the full field
    p.taper=zeros(p.NyNx);
    qy=floor(p.NyNx(1)/4);qx=floor(p.NyNx(2)/4);
    p.taper(qy:3*qy-1,qx:3*qx-1)=deal(1.);
    p.mask=p.taper; 
    [v,I,w,vw,scl,cr]=maskit(Hx,p,scl,Hx);

    % Number of iterations
    defval('xver',100);

    %limspc=0.35;
    %datspc=0.1;

    % Labels for the nu axes will be based on the special values
    nutickvals=[1/3 1/2 1 3/2 5/2];
    %nuticklabs={'$1/3$','$1/2$','$1$','$\nu_0$','$5/2$'};
    nuticklabs={'$\frac{1}{3}$','$\frac{1}{2}$','$1$','$\nu_0=\frac{3}{2}$','$\frac{5}{2}$'};
    % L=floor(sqrt(prod(p.NyNx)));
    L=floor(sqrt(sum((p.NyNx.*p.dydx).^2)));
    % color options for kde plots
    %clrs=["#D95F02","#0072BD","#EDB120","#7570B3","#1B9E77"];
    clrs=[[0.6 0.6 0.6],"7570B3","#0072BD","#EDB120","#1B9E77"];
    % font size
    afs=10;
    % marker size
    ms=100;
    % power of 2 for kde
    kdep=8;
    
    defval('fixflag',0)
    thftrial=[th(1) 5/2 th(3)];

    % The full field
    %filesuff=append('fullrect_',datum);
    filesuff=append('fullrect_',datum);
    try
        thhat_r1=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r1,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_r1,[],filesuff);
        thhat_r1=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r1,1)<xver 
            mleosl('demo1',xver,th,p_r1,[],filesuff);
            thhat_r1=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    if fixflag==0
        filesuff=append(filesuff,'_f');
        thf=th;
    elseif fixflag==1
        filesuff=append(filesuff,'_f_INCORR');
        thf=thftrial;
    end
    try
        thhat_r1_f=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r1_f,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_r1,[],filesuff,thf,[1 0 1]);
        thhat_r1_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r1_f,1)<xver
            mleosl('demo1',xver,th,p_r1,[],filesuff,thf,[1 0 1]);
            thhat_r1_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The padded rectangle case
    p_r2=p;
    filesuff=append('padrect_',datum);
    try
        thhat_r2=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r2,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_r2,[],filesuff);
        thhat_r2=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r2,1)<xver
            mleosl('demo1',xver,th,p_r2,[],filesuff);
            thhat_r2=load(append('mleosl_thhat_',filesuff));
        end
    end
    % including for fixed nu
    if fixflag==0
        filesuff=append(filesuff,'_f');
        thf=th;
    elseif fixflag==1
        filesuff=append(filesuff,'_f_INCORR');
        thf=thftrial;
    end
    try
        thhat_r2_f=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r2_f,1)<xver;error('Calculating additional realizations');end
    catch
        %filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th,p_r2,[],filesuff,thf,[1 0 1]);
        thhat_r2_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r2_f,1)<xver 
            mleosl('demo1',xver,th,p_r2,[],filesuff,thf,[1 0 1]);
            thhat_r2_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The small rectangle case
    p_r3=p;
    qy=floor(p_r3.NyNx(1)/4);qx=floor(p_r3.NyNx(2)/4);
    p_r3.NyNx=[2*qy 2*qx];
    p_r3.taper=ones(p_r3.NyNx);
    p_r3.mask=p_r3.taper; 
    filesuff=append('smallrect_',datum);
    try
        thhat_r3=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r3,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th,p_r3,[],filesuff);
        thhat_r3=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r3,1)<xver
            mleosl('demo1',xver,th,p_r3,[],filesuff);
            thhat_r3=load(append('mleosl_thhat_',filesuff));
        end
    end

    % including for fixed nu
    if fixflag==0
        filesuff=append(filesuff,'_f');
        thf=th;
    elseif fixflag==1
        filesuff=append(filesuff,'_f_INCORR');
        thf=thftrial;
    end
    try
        thhat_r3_f=load(append('mleosl_thhat_',filesuff));
        if size(thhat_r3_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff);
        mleosl('demo1',xver,th,p_r3,[],filesuff,thf,[1 0 1]);
        thhat_r3_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_r3_f,1)<xver 
            mleosl('demo1',xver,th,p_r3,[],filesuff,thf,[1 0 1]);
            thhat_r3_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    for hnd=1:3
       if hnd==1
       % The full field -- we expect the estimates to have less variance,
       % narrower histograms
           fld=vw;mt=0;
           thhat1=thhat_r1;thhat1f=thhat_r1_f;
       elseif hnd==2
       % The inner rectangle -- we expect that the padding by zeros case shoule
       % be equivalent to just having subfield
           fld=zeros(p.NyNx);vw2=v2s(vw,p);
           fld(qy:3*qy-1,qx:3*qx-1)=vw2(qy:3*qy-1,qx:3*qx-1);fld=fld(:);mt=1;
           thhat1=thhat_r2;thhat1f=thhat_r2_f;
       elseif hnd==3
       % The subfield
           vw2=v2s(vw,p);fld=vw2(qy:3*qy-1,qx:3*qx-1);p.NyNx=size(fld);mt=0;
           thhat1=thhat_r3;thhat1f=thhat_r3_f;
       end
        % The estimates associated with the current column
        data1=thhat1(1:xver,:);data1f=thhat1f(1:xver,:);
        t1m=mean(data1);t1fm=mean(data1f);
        t1sd(hnd,:)=std(data1);t1fsd(hnd,:)=std(data1f);

        % The random field plot
        axes(ha((hnd-1)*numrows+1))
        [h,t]=plotit(fld,p,cr,[],mt);
        longticks
        if hnd==2
            set(h.Parent,'YTickLabel',[]);
            b=rectangle('Position',[qx, qy, 2*qx-1, 2*qy-1],...
                'EdgeColor','k','LineWidth',2);
        elseif hnd==3
            shrink(ha((hnd-1)*numrows+1),2,2);
            b=rectangle('Position',[1, 1, p.NyNx(2)-1, p.NyNx(1)-1],...
                'EdgeColor','k','LineWidth',2);
        end
        set(gca,'FontSize',afs,'XTickLabelRotationMode','manual');

        % The variance estimates
        axes(ha((hnd-1)*numrows+2))
        hold on
        datspc=th(1)/2;
        kdelim=[min(data1(:,1)) max(data1(:,1))]+[-1 1].*datspc;
        [a11,bdens11,c11]=kde(data1(:,1),2^kdep,kdelim(1),kdelim(2));
        %t11kp=bar(c11,t1sd(hnd,1)*bdens11,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        t11kp=bar(c11,bdens11,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        kdelim=[min(data1f(:,1)) max(data1f(:,1))]+[-1 1].*datspc;
        [a11f,bdens11f,c11f]=kde(data1f(:,1),2^kdep,kdelim(1),kdelim(2));
        %t11fkp=bar(c11f,-1*t1fsd(hnd,1)*bdens11f,1,'FaceColor',clrs{2},'EdgeColor',clrs{2});
        t11fkp=bar(c11f,-1*bdens11f,1,'FaceColor',clrs{2},'EdgeColor',clrs{2});
        th01=scatter(th(1),0,ms,'k|');
        %maxy(hnd,1)=max(abs([t1sd(hnd,1)*bdens11; t1fsd(hnd,1)*bdens11f]));
        maxy(hnd,1)=max(abs([bdens11; bdens11f]));
        maxx(hnd,1)=max(abs([c11 c11f]));
        plot([min(abs([c11 c11f])),maxx(hnd,1)],[0,0],'k')
        %xlim(th(1)+limspc.*max(abs([data1(:,1); data1f(:,1)]))*[-1 1])
        %ylim([-1*miny maxy].*1.3)
        if hnd==1
            ylabel('$\hat{\sigma}^2$','Interpreter', 'latex')
            set(get(gca,'ylabel'),'rotation',0)
        end
        longticks
        set(gca,'YTickLabel',[])
        % xtickvals=gca().XAxis.TickValues;
        % xtickvals=[xtickvals(1) xtickvals(end)];
        % if any(xtickvals==th(1));
        %     sigticklabs=cellstr(num2str(xtickvals(:)./L^2,'%0.2g'));
        %     sigticklabs{xtickvals==th(1)}=sprintf('$\\frac{\\sigma^2_0}{[L^2]}$=%0.2g',th(1)/L^2);
        % else
        %     [xtickvals,srtsig]=sort([xtickvals th(1)]);
        %     sigticklabs=cellstr(num2str(xtickvals(:)./L^2,'%0.2g'));
        %     sigticklabs{srtsig==max(srtsig)}=sprintf('$\\frac{\\sigma^2_0}{[L^2]}$=%0.2g',th(1)/L^2);
        % end
        % set(gca,'XTick',xtickvals,'XTickLabel',sigticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');
        box on
        hold off

        % The smoothness estimates
        axes(ha((hnd-1)*numrows+3))
        hold on
        datspc=th(2)/2;
        kdelim=[min(data1(:,2)) max(data1(:,2))]+[-1 1].*datspc;
        [a12,bdens12,c12]=kde(data1(:,2),2^kdep,kdelim(1),kdelim(2));
        %t12kp=bar(c12,t1sd(hnd,2)*bdens12,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        t12kp=bar(c12,bdens12,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        th02=scatter(th(2),0,ms,'k|');
        th02f=scatter(thf(2),-1.0,floor(ms/2),'kv');
        %maxy(hnd,2)=max(abs([t1sd(hnd,2)*bdens12]));
        maxy(hnd,2)=max(abs([bdens12]));
        maxx(hnd,2)=max(abs(c12));
        plot([min(abs(c12)),maxx(hnd,2)],[0,0],'k')
        if hnd==1
            ylabel('$\hat{\nu}$','Interpreter', 'latex')
            set(get(gca,'ylabel'),'rotation',0)
        end
        longticks
        set(gca,'YTickLabel',[])
        xtickvals=gca().XAxis.TickValues;
        set(gca,'XTick',nutickvals,'XTickLabel',nuticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');
        box on
        hold off

        % The range estimates
        axes(ha((hnd-1)*numrows+4))
        hold on
        datspc=th(3)/2;
        kdelim=[min(data1(:,3)) max(data1(:,3))]+[-1 1].*datspc;
        %[a13,bdens13,c13]=kde(data1(:,3),2^kdep);
        [a13,bdens13,c13]=kde(data1(:,3),2^kdep,kdelim(1),kdelim(2));
        %t13kp=bar(c13,t1sd(hnd,3)*bdens13,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        t13kp=bar(c13,bdens13,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        kdelim=[min(data1f(:,3)) max(data1f(:,3))]+[-1 1].*datspc;
        %[a13f,bdens13f,c13f]=kde(data1f(:,3),2^kdep);
        [a13f,bdens13f,c13f]=kde(data1f(:,3),2^kdep,kdelim(1),kdelim(2));
        %t13fkp=bar(c13f,-1*t1fsd(hnd,3)*bdens13f,1,'FaceColor',clrs{2},'EdgeColor',clrs{2});
        t13fkp=bar(c13f,-1*bdens13f,1,'FaceColor',clrs{2},'EdgeColor',clrs{2});
        th03=scatter(th(3),0,ms,'k|');
        %maxy(hnd,3)=max(abs([t1sd(hnd,3)*bdens13; t1fsd(hnd,3)*bdens13f]));
        maxy(hnd,3)=max(abs([bdens13; bdens13f]));
        maxx(hnd,3)=max(abs([c13 c13f]));
        plot([min(abs(c13)),maxx(hnd,3)],[0,0],'k')
        if hnd==1
            ylabel('$\hat{\rho}$','Interpreter', 'latex')
            set(get(gca,'ylabel'),'rotation',0)
        end
        longticks
        set(gca,'YTickLabel',[])
        box on
        hold off
    end

    % Now assign a common axis range for all of the kde subplots
    %xlspc=1;
    ylspc=1.1;
    stdOut=5;
    for hnd=1:3
    %    set(H(hnd+1,:),'XLim',th(hnd)+(max(maxx(:,hnd))-th(hnd)).*[-1*xlspc xlspc],...
         set(H(hnd+1,:),'XLim',th(hnd)+max(t1sd(:,hnd)).*[-stdOut stdOut],...
          'YLim',max(maxy(:,hnd)).*[-1*ylspc ylspc]);
    end
    % Treat labels for sigma and rho
    for hnd=1:3
        axes(H(2,hnd));
        xtickvals=gca().XAxis.TickValues;
        if any(xtickvals==th(1));
            xticklabs=cellstr(num2str(xtickvals(:),'%0.1g'));
            xticklabs{xtickvals==th(1)}=sprintf('$\\sigma^2_0$=%0.1g',th(1));
            xticklabs={xticklabs{1}; xticklabs{xtickvals==th(1)};...
                xticklabs{end}};
            xtickvals=[xtickvals(1) th(1) xtickvals(end)];
        else
            [xtickvals,srtsig]=sort([xtickvals th(1)]);
            xticklabs=cellstr(num2str(xtickvals(:),'%0.1g'));
            xticklabs{srtsig==max(srtsig)}=sprintf('$\\sigma^2_0$=%0.1g',th(1));
            xticklabs={xticklabs{1}; xticklabs{srtsig==max(srtsig)};...
                xticklabs{end}};
            xtickvals=[xtickvals(1) th(1) xtickvals(end)];
        end
        set(gca,'XTick',xtickvals,'XTickLabel',xticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');

        axes(H(4,hnd));
        xtickvals=gca().XAxis.TickValues;
        if any(xtickvals==th(3));
            xticklabs=cellstr(num2str(xtickvals(:)./L,'%0.1g'));
            xticklabs{xtickvals==th(3)}=sprintf('$\\rho_0/[L]$=%0.1g',th(3)/L);
            xticklabs={xticklabs{1}; xticklabs{xtickvals==th(3)};...
                xticklabs{end}};
            xtickvals=[xtickvals(1) th(3) xtickvals(end)];
        else
            [xtickvals,srtrh]=sort([xtickvals th(3)]);
            xticklabs=cellstr(num2str(xtickvals(:)./L,'%0.1g'));
            xticklabs{srtrh==max(srtrh)}=sprintf('$\\rho_0/[L]$=%0.1g',th(3)/L);
            xticklabs={xticklabs{1}; xticklabs{srtrh==max(srtrh)};...
                xticklabs{end}};
            xtickvals=[xtickvals(1) th(3) xtickvals(end)];
        end
        set(gca,'XTick',xtickvals,'XTickLabel',xticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');
    end
    % Special treatment for nu x-axis
    xlspc=1.1;
    set(H(3,:),'XLim',th(2)+(max(abs(th(2)-nutickvals))).*[-1*xlspc xlspc]);
    %[0.1 2.73]

    % Shift
    axes(H(1,1));rf1pos=gca().Position;axes(H(1,2));rf2pos=gca().Position;
    axes(H(1,3));rf3pos=gca().Position;
    fshr=0.65;
    rf1npos=shrink(H(1,1),fshr,fshr);
    rf2npos=shrink(H(1,2),fshr,fshr);
    rf3npos=shrink(H(1,3),fshr,fshr);
    H(1,1).Position(2)=rf1pos(2)-rf1npos(4)+rf1pos(4);
    H(1,2).Position(2)=rf2pos(2)-rf2npos(4)+rf2pos(4);
    H(1,3).Position(2)=rf3pos(2)-rf3npos(4)+rf3pos(4);

    thshr=1.5;
    shrink([H(2:4,:)],1,thshr)
    serre([H(2:4,1)],1/2,'down');serre([H(2:4,2)],1/2,'down');serre([H(2:4,3)],1/2,'down')
    serre([H(1,:)' H(2,:)'],1/4,'down');

    % Minimize space between the subplots
    serre([H(2:4,:)],1/2,'across')
    serre([H(1,:)],3/4,'across');serre([H(1,:)],3/4,'across');serre([H(1,:)],3/4,'across'); 
    % Save once we are happy
    fname=sprintf('samplingfigA1_%3.3i_%3.3i_%i_fld',p.NyNx,xver);
    if fixflag==1;fname=append(fname,'_INCORR');end;
    figdisp([],fname,[],1);
end

function varargout=plotit(v,p,cr,th,mt,afs)
% use for consistent imagefnan calls NaN fill color:
% nancol=abyss(2);
% h=imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],halverange(v,80,NaN),nancol(1,:)); 
h=imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),'parula',halverange(v,80,NaN));
axis ij image
%hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
defval('mt',1)
defval('afs',10);
if mt==1
  if ~isempty(th)
    if isstr(th)
        t=title(th);
    else
        t=title(sprintf('%i x %i | %s = [%g %g %g]',p.NyNx(1),p.NyNx(2),'\theta',...
                        th./[1 1 sqrt(prod(p.dydx))]));
    end
    t.FontSize=afs+1;
    movev(t,-p.NyNx(1)/25)
  else
    t=[];
  end
else
  t=[];
end
hold on
plot([1 p.NyNx(2)],[1 1],'k')
plot([1 p.NyNx(2)],[p.NyNx(1) p.NyNx(1)],'k')
plot([1 1],[1 p.NyNx(1)],'k')
plot([p.NyNx(2) p.NyNx(2)],[1 p.NyNx(1)],'k')
varns={h,t};
varargout=varns;
