function specialcases(th,params,est,xver,fixflag)
% SPECIALCASES(TH,P,EST)
%
% Calculates and renders paper figures for demonstrating special smoothness cases
%
% Input:
%
% th        Matern parameters
% params    Grid parameters
% est       Couple versions of plot possible, 0 and 2 are current
% xver      Number of iterations, borrowing naming convention from BLUROSY
% fixflag   Flag to set the fixed case to the WRONG nu, 5/2; otherwise, fix to
%           correct nu; only relevant when est==2
%
% Output:
%
% <>        Generates and saves a figure (and possibly runs and stores many estimate
%           calculations)
% 
% Examples:
%
% Figure 1: Illustration of the special cases
%
% specialcases([],[],0)
%
% For estimation with nu fixed to the wrong value:
% specialcases([],[],2,[],1)
%
% For estimation with nu fixed to the correct value (NEED TO RUN THIS):
% specialcases([],[],2,[],0)
%
% Last modified by olwalbert-at-princeton.edu on 02/09/2026

% Flag for Figure 4 or Figure 5
defval('est',0)
% Field size 
defval('params',[188 233]+randi(20,[1 2]))
% Rely on default scaling defined within MUCKIT or MASKIT
defval('scl',[]);
% Some combinations - SIMULOSL
%defval('th',[1 1.5 5e3]);
defval('th',[1 1.5 3]);
p.NyNx=[params(1) params(2)];
p.dydx=[1 1];
%p.dydx=1e3*[1 1];
p.blurs=Inf;
p.blurs=0;

%datum=date;
datum='05-Sep-2024'; %use for final figures -- need to figure out how to
%                    display the incorrectly fixed values better. for
%                    demonstration purposes, it might be worth choosing a
%                    smoothness of 1, closer to 1/3 and 1/2. Infinity case
%                    should receive its own axis range 
%datum='21-Jun-2024'; %used for debugging; there were errors in fixing nu for estimation in this
%                     suite

if est==0
    % Create a figure in the style of Fig 4 of the sampling and special cases
    % paper, demonstrating the 6 special cases as simulations only
    p.NyNx=[122,141];
    th1=[1 1/3 5];
    th2=[1 1/2 5];
    th3=[1 1 5];
    th4=[1 3/2 5];
    th5=[1 5/2 5];
    th6=[1 Inf 5];

    labs={'$$\nu = 1/{}3$$',%\frac{1}{3}$$',
          '$$\nu = 3/{}2$$',%\frac{3}{2}$$',
          '$$\nu = 1/{}2$$',%\frac{1}{2}$$',
          '$$\nu = 5/{}2$$',%\frac{5}{2}$$',
          '$$\nu = 1$$',
          '$$\nu \rightarrow \infty$$'};

    figure(1);clf
    set(gcf,'InnerPosition',[1 1 1140 760])
    [ah,ha,H]=krijetem(subnum(2,3));
    for hnd=1:3
        if hnd==1
            % The one-third and three-halves case 
            [Hx,~,~,~,~,~,~,~,~,Z1]=simulosl(th1,p);nu=th1(2);
            Hx2=simulosl(th4,p,[],Z1);nu2=th4(2);
            mt=0;
        elseif hnd==2
            % The one-half and five-halves case 
            Hx=simulosl(th2,p,[],Z1);nu=th2(2);
            Hx2=simulosl(th5,p,[],Z1);nu2=th5(2);
            mt=0;
        elseif hnd==3
            % The 1 and Inf case 
            Hx=simulosl(th3,p,[],Z1);nu=th3(2);
            Hx2=simulosl(th6,p,[],Z1);nu2=th6(2);
            mt=0;
            %mt=1;
        else                                                                    
            error('there is an indexing error; debug creation of subfigures')   
        end                                                                     
        axes(ha((hnd-1)*2+1))                                                   
        [h,cax,t,cb]=plotit(Hx,p,[],th,mt,ha((hnd-1)*2+1));
        axes(ha((hnd-1)*2+1))                                                   
        st(hnd)=title(sprintf('%s',labs{(hnd-1)*2+1}),'interpreter', 'latex');
        longticks
        set(h.Parent,'XTicklabel',[]);
        if hnd>1
            set(h.Parent,'YTickLabel',[]);
        end
        axes(ha((hnd-1)*2+2))                                                   
        [h2,cax2,t2,cb2]=plotit(Hx2,p,[],th,0,ha((hnd-1)*2+2));
        axes(ha((hnd-1)*2+2))                                                   
        st2(hnd)=title(sprintf('%s',labs{(hnd-1)*2+2}),'interpreter', 'latex');
        longticks
        if hnd>1
            set(h2.Parent,'YTickLabel',[]);
        end
        if hnd<3
            delete([cb cb2])
        end
    end  
    % all of the objects that must be moved
    ho2=getpos(ha);o2=serre(H',1,'down');
    movev(H(1,:),-0.01);
    o3=serre(H,1/2,'across');
    cbp2=cb.Position(2);
    cb.Position(2)=cbp2-(-o3(1,2)+o2(1,2));
    cb2p2=cb2.Position(2);
    cb2.Position(2)=cb2p2+(-o3(1,2)+o2(1,2))/4;
    movev(H,0.05);moveh(H,-0.05);
    movev([cb cb2],0.05);moveh([cb cb2],-0.05);
    movev([st st2],4)
    % some more fixing
    moveh(cb2,0.01)
    moveh(cb,0.01)
    for ind=1:6
      ah(ind).Title.FontSize=12.5;
    end
    ah(1).YLabel.String='northing [m]'
    ah(4).YLabel.String='northing [m]';
    ah(4).XLabel.String='easting [m]';
    ah(5).XLabel.String='easting [m]';
    ah(6).XLabel.String='easting [m]';
    keyboard
    figdisp([],sprintf('specialfig4_%3.3i_%3.3i_fld',p.NyNx),[],1);
elseif est==1
    % Create a figure in the style of Fig 5 of the sampling and special cases
    % paper, demonstrating common special cases as simulations, their estimates
    % with uncertainty for 3 and 2 parameter cases, where the 2 parameter case
    % is fixed to the WRONG nu, 5/2
    
    % We will show a panel for (1) a von karman (nu=1/3) field, 
    % (2) an exponential (nu=1/2) field, and (3) a Gaussian (nu->Inf) field
    % Number of iterations
    defval('xver',100);

    % Once we find a grid-size we like through experimentation, let's set this
    % size for future program calls
    p.NyNx = [188, 233]; 
    th1 = [1 1/3 5];
    th2 = [1 1/2 5];
    th3 = [1 Inf 5];
    th4 = [1 5/2 5];
    
    labs={'$$\nu_0 = \frac{1}{3}$$',
          '$$\nu_0 = \frac{1}{2}$$',
          '$$\nu_0 \rightarrow \infty$$'};

    % Nu=1/3 case
    filesuff=append('nuthird_',datum);
    try
        thhat_1=load(append('mleosl_thhat_',filesuff));
        if size(thhat_1,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th1,p,[],filesuff);
        thhat_1=load(append('mleosl_thhat_',filesuff));
        while size(thhat_1,1)<xver 
            mleosl('demo1',xver,th1,p,[],filesuff);
            thhat_1=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu=1/3 case fixed to nu=5/2
    try
        thhat_1_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_1_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th1,p,[],filesuff,th4,[1 0 1]);
        thhat_1_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_1_f,1)<xver 
            mleosl('demo1',xver,th1,p,[],filesuff,th4,[1 0 1]);
            thhat_1_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu=1/2 case
    filesuff=append('nuhalf_',datum);
    try
        thhat_2=load(append('mleosl_thhat_',filesuff));
        if size(thhat_2,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th2,p,[],filesuff);
        thhat_2=load(append('mleosl_thhat_',filesuff));
        while size(thhat_2,1)<xver 
            mleosl('demo1',xver,th2,p,[],filesuff);
            thhat_2=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu=1/2 case fixed to nu=5/2
    try
        thhat_2_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_2_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th2,p,[],filesuff,th4,[1 0 1]);
        thhat_2_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_2_f,1)<xver 
            mleosl('demo1',xver,th2,p,[],filesuff,th4,[1 0 1]);
            thhat_2_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu->Inf case; note that we want to estimate nu freely in the 3 parameter
    % case, this means that we cannot set thini or aguess 
    filesuff=append('nuinf_',datum);
    try
        thhat_3=load(append('mleosl_thhat_',filesuff));
        if size(thhat_3,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th3,p,[],filesuff);
        thhat_3=load(append('mleosl_thhat_',filesuff));
        while size(thhat_3,1)<xver 
            mleosl('demo1',xver,th3,p,[],filesuff);
            thhat_3=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu->Inf case fixed to nu=5/2
    try
        thhat_3_f=load(append('mleosl_thhat_',filesuff,'_f'));
        if size(thhat_3_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f');
        mleosl('demo1',xver,th3,p,[],filesuff,th4,[1 0 1]);
        thhat_3_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_3_f,1)<xver 
            mleosl('demo1',xver,th3,p,[],filesuff,th4,[1 0 1]);
            thhat_3_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The setup for creating the plots of the simulated fields is very similar
    % to est==0, but we will need to load our previously stored  XVER estimates
    % from their MLEOSL hash
    
    % maxsd needs to be calculated for every set of thhats:  
    maxsd=max(max([std(thhat_1);std(thhat_1_f);std(thhat_2);std(thhat_2_f);
               std(thhat_3);std(thhat_3_f)]));
    p2=p;
    symbs={'^','s','o'};
    numrows=3;
    figure(1)
    clf
    [ah,ha,H]=krijetem(subnum(3,3));
    for hnd=1:3
        if hnd==1
            % The one-third case 
            Hx=simulosl(th1,p);
            thhatf=NaN;
            while isnan(thhatf);
                [thhatf,~,~,sclf]=mleosl(Hx,[],p,[],[],th4,[1 0 1],0);
            end
            thhatf=thhatf.*sclf
            Hx2=simulosl(thhatf,p);
            thhat1=thhat_1;thhat1f=thhat_1_f;
            mt=0;
        elseif hnd==2
            % The one-half case 
            Hx=simulosl(th2,p);
            thhatf=NaN;
            while isnan(thhatf);
                [thhatf,~,~,sclf]=mleosl(Hx,[],p,[],[],th4,[1 0 1],0);
            end
            thhatf=thhatf.*sclf
            Hx2=simulosl(thhatf,p);
            thhat1=thhat_2;thhat1f=thhat_2_f;
            mt=0;
        elseif hnd==3
            % The Inf case 
            Hx=simulosl(th3,p);
            thhatf=NaN;
            while isnan(thhatf);
                [thhatf,~,~,sclf]=mleosl(Hx,[],p,[],[],th4,[1 0 1],0);
            end
            try
                Hx2=simulosl(thhatf.*sclf,p);
            catch
                p2=p;p2.NyNx=p2.NyNx.*4;
                try                                                                 
                    Hx2_=simulosl(thhatf.*sclf,p2);                                   
                    Hx2=quarter(v2s(Hx2_,p2));
                catch    
                    p2.NyNx=p2.NyNx.*4;
                    try                                                             
                        Hx2_=simulosl(thhatf.*sclf,p2);                              
                        Hx2=quarter(quarter(v2s(Hx2_,p2)));
                    catch                                                           
                        p2.NyNx=p2.NyNx.*4; 
                        Hx2_=simulosl(thhatf.*sclf,p2);
                        Hx2=quarter(quarter(quarter(v2s(Hx2_,p2))));
                    end
                end
            end
            thhatf=thhatf.*sclf;
            thhat1=thhat_3;thhat1f=thhat_3_f;
            mt=0;
        else                                                                    
            error('there is an indexing error; debug creation of subfigures')   
        end                                                                     
        axes(ha((hnd-1)*numrows+1))                                                   
        t=plotit(Hx,p,[],th,mt);
        if ~isempty(t)
            t.Position(1)=-3*t.Position(1);
            %t.Position(1)=-1.75*t.Position(1);
        end
        st=title(sprintf('%s',labs{hnd}),'interpreter', 'latex');

        axes(ha((hnd-1)*numrows+2))
        % Scatter marker for each th0 value, scale will be standardized by +/-3
        % std dev
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
        thhat1fsp=plot([t1fm-t1fsd;t1fm+t1fsd],[yfvls(1,:);yfvls(1,:)],'r-');
        thhat1sp=plot([t1m-t1sd;t1m+t1sd],[yvls(1,:);yvls(1,:)],'k-');
        thhat1fspb=plot([t1fm-t1fsd;t1fm+t1fsd],[yfvls(1,:);yfvls(1,:)],'r|');
        thhat1spb=plot([t1m-t1sd;t1m+t1sd],[yvls(1,:);yvls(1,:)],'k|');
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

        axes(ha((hnd-1)*numrows+3))                                                   
        plotit(Hx2,p2,[],th,0);
    end  
    figdisp([],sprintf('specialfig5_%3.3i_%3.3i_%3.3i_fld',p.NyNx,xver),[],1);
elseif est==2
    % Create a figure in the style of Fig 5 of the sampling and special cases
    % paper, demonstrating common special cases as simulations, their estimates
    % with uncertainty for 3 and 2 parameter cases, where the 2 parameter case
    % is fixed, either to the correct (fixflag=0) or to the WRONG 
    % (fixflag=1) nu, 5/2; this figure has improved plots of parameter
    % estimates shown in a violin style
    defval('fixflag',0);
    % We will show a panel for (1) a von karman (nu=1/3) field, 
    % (2) an exponential (nu=1/2) field, and (3) a Gaussian (nu->Inf) field
    % Number of iterations
    defval('xver',200);

    % Once we find a grid-size we like through experimentation, let's set this
    % size for future program calls
    p.NyNx = [188, 233]; 
    th1 = [1 1/3 3];
    th2 = [1 1/2 3];
    th3 = [1 Inf 3];
    th4 = [1 3/2 3]; % previously, nu=5/2;
    
    labs={sprintf('$$\\theta_0 = $$[%g 1/3 %g]',th1(1),th1(3)),
          sprintf('$$\\theta_0 = $$[%g 1/2 %g]',th2(1),th2(3)),
          sprintf('$$\\theta_0 = $$[%g $$\\infty$$ %g]',th3(1),th3(3))};

    % Nu=1/3 case
    filesuff=append('nuthird_',datum);
    try
        thhat_1=load(append('mleosl_thhat_',filesuff));
        if size(thhat_1,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th1,p,[],filesuff);
        thhat_1=load(append('mleosl_thhat_',filesuff));
        while size(thhat_1,1)<xver 
            mleosl('demo1',xver,th1,p,[],filesuff);
            thhat_1=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu=1/3 case fixed
    try
        thhat_1_f=load(append('mleosl_thhat_',filesuff,'_f',num2str(fixflag)));
        if size(thhat_1_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f',num2str(fixflag));
        if fixflag==0;thf=th1;elseif fixflag==1;thf=th4;end
        mleosl('demo1',xver,th1,p,[],filesuff,thf,[1 0 1]);
        thhat_1_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_1_f,1)<xver 
            mleosl('demo1',xver,th1,p,[],filesuff,thf,[1 0 1]);
            thhat_1_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu=1/2 case
    filesuff=append('nuhalf_',datum);
    try
        thhat_2=load(append('mleosl_thhat_',filesuff));
        if size(thhat_2,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th2,p,[],filesuff);
        thhat_2=load(append('mleosl_thhat_',filesuff));
        while size(thhat_2,1)<xver 
            mleosl('demo1',xver,th2,p,[],filesuff);
            thhat_2=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu=1/2 case fixed
    try
        thhat_2_f=load(append('mleosl_thhat_',filesuff,'_f',num2str(fixflag)));
        if size(thhat_2_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f',num2str(fixflag));
        if fixflag==0;thf=th2;elseif fixflag==1;thf=th4;end
        mleosl('demo1',xver,th2,p,[],filesuff,thf,[1 0 1]);
        thhat_2_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_2_f,1)<xver 
            mleosl('demo1',xver,th2,p,[],filesuff,thf,[1 0 1]);
            thhat_2_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu->Inf case; note that we want to estimate nu freely in the 3 parameter
    % case, this means that we cannot set thini or aguess 
    filesuff=append('nuinf_',datum);
    try
        thhat_3=load(append('mleosl_thhat_',filesuff));
        if size(thhat_3,1)<xver;error('Calculating additional realizations');end
    catch
        mleosl('demo1',xver,th3,p,[],filesuff);
        thhat_3=load(append('mleosl_thhat_',filesuff));
        while size(thhat_3,1)<xver 
            mleosl('demo1',xver,th3,p,[],filesuff);
            thhat_3=load(append('mleosl_thhat_',filesuff));
        end
    end

    % Nu->Inf case fixed
    try
        thhat_3_f=load(append('mleosl_thhat_',filesuff,'_f',num2str(fixflag)));
        if size(thhat_3_f,1)<xver;error('Calculating additional realizations');end
    catch
        filesuff=append(filesuff,'_f',num2str(fixflag));
        if fixflag==0;thf=th3;elseif fixflag==1;thf=th4;end
        mleosl('demo1',xver,th3,p,[],filesuff,thf,[1 0 1]);
        thhat_3_f=load(append('mleosl_thhat_',filesuff));
        while size(thhat_3_f,1)<xver 
            mleosl('demo1',xver,th3,p,[],filesuff,thf,[1 0 1]);
            thhat_3_f=load(append('mleosl_thhat_',filesuff));
        end
    end

    % The setup for creating the plots of the simulated fields is very similar
    % to est==0, but we will need to load our previously stored  XVER estimates
    % from their MLEOSL hash
    
    % maxsd needs to be calculated for every set of thhats:  
    maxsd=max(max([std(thhat_1);std(thhat_1_f);std(thhat_2);std(thhat_2_f);
               std(thhat_3);std(thhat_3_f)]));
    p2=p;
    symbs={'^','s','o'};
    numrows=5;
    datspc=0.1;
    limspc=0.35;
    L=floor(sqrt(prod(p.NyNx)));
    % color options for kde plots
    clrs=["#D95F02","#0072BD","#EDB120","#7570B3","#1B9E77"];
    % font size
    afs=10;
    % marker size
    ms=100;
    % power of 2 for kde
    kdep=8;
    figure(4)
    clf
    [ah,ha,H]=krijetem(subnum(numrows,3));
    for hnd=1:3
        if hnd==1
            % The one-third case 
            th=th1;
            if fixflag==0;thf=th1;elseif fixflag==1;thf=th4;end
            Hx=simulosl(th1,p);
            dbgolw=0;
            if dbgolw==1
                thhatf=th1;  
                Hx2=Hx;
            else
                thhatf=median(thhat_1)
                %thhat_1(randi(size(thhat_1,1)),:);
                %thhatf=NaN;
                %while isnan(thhatf);
                %    [thhatf,~,~,sclf]=mleosl(Hx,[],p,[],[],thf,[1 0 1],0);
                %end
                %thhatf=thhatf.*sclf
                Hx2=simulosl(thhatf,p);
            end
            thhat1=thhat_1;thhat1f=thhat_1_f;
            mt=0;
        elseif hnd==2
            % The one-half case 
            th=th2;
            if fixflag==0;thf=th2;elseif fixflag==1;thf=th4;end
            Hx=simulosl(th2,p);
            if dbgolw==1
                thhatf=th1;  
                Hx2=Hx;
            else
                thhatf=median(thhat_2)
                %thhatf=thhat_2(randi(size(thhat_2,1)),:);
                Hx2=simulosl(thhatf,p);
            end
            thhat1=thhat_2;thhat1f=thhat_2_f;
            mt=0;
        elseif hnd==3
            th=th3;
            if fixflag==0;thf=th3;elseif fixflag==1;thf=th4;end
            % The Inf case 
            Hx=simulosl(th3,p);
            if dbgolw==1
                thhatf=th1;  
                Hx2=Hx;
            else
                thhatf=median(thhat_3)
                %thhatf=thhat_3(randi(size(thhat_3,1)),:);
                try
                    Hx2=simulosl(thhatf,p);
                catch
                    p2=p;p2.NyNx=p2.NyNx.*4;
                    try                                                       
                        Hx2_=simulosl(thhatf,p2);
                        Hx2=Hx2_;
                        %Hx2=quarter(quarter(v2s(Hx2_,p2)));
                        %Hx2=Hx2(:);p2=p;
                    catch    
                        p2.NyNx=p2.NyNx.*4;
                        try                                                             
                            Hx2_=simulosl(thhatf,p2);
                            Hx2=Hx2_;
                            %Hx2=quarter(quarter(quarter(quarter(v2s(Hx2_,p2)))));
                            %Hx2=Hx2(:);p2=p;
                        catch                                                           
                            p2.NyNx=p2.NyNx.*4; 
                            Hx2_=simulosl(thhatf,p2);
                            Hx2=Hx2_;
                            %Hx2=quarter(quarter(quarter(quarter(quarter(quarter(v2s(Hx2_,p2)))))));
                            %Hx2=Hx2(:);p2=p;
                        end
                    end
                end
            end
            thhat1=thhat_3;thhat1f=thhat_3_f;
            mt=0;
        else                                                                    
            error('there is an indexing error; debug creation of subfigures')   
        end                                                                     
        % standardize estimates by centering on the mean and scaling the kde by
        % the observed standard deviation::
        data1=thhat1(1:xver,:);data1f=thhat1f(1:xver,:);
        t1m=mean(data1);t1fm=mean(data1f);
        t1sd=std(data1);t1fsd=std(data1f);

        axes(ha((hnd-1)*numrows+1))
        hold on
        [h,cax,t]=plotit(Hx,p,[],th,mt);
        longticks
        if hnd>1;set(gca,'YTickLabel',[]);end
        set(gca,'XTickLabel',[],'FontSize',9);
        st(hnd)=title(sprintf('%s',labs{hnd}),'interpreter', 'latex','FontSize',afs);
        hold off

        % The variance estimates
        axes(ha((hnd-1)*numrows+2))
        hold on
        datspc=th(1)/2;
        kdelim=[min(data1(:,1)) max(data1(:,1))]+[-1 1].*datspc;
        [a11,bdens11,c11]=kde(data1(:,1),2^kdep,kdelim(1),kdelim(2));
        t11kp=bar(c11,bdens11,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        kdelim=[min(data1f(:,1)) max(data1f(:,1))]+[-1 1].*datspc;
        [a11f,bdens11f,c11f]=kde(data1f(:,1),2^kdep,kdelim(1),kdelim(2));
        t11fkp=bar(c11f,-1*bdens11f,1,'FaceColor',clrs{2},'EdgeColor',clrs{2});
        th01=scatter(th(1),0,ms,'k|');
        thf1=scatter(thhatf(1),0,floor(ms/4),'ko');
        maxx(hnd,1)=max(abs([c11 c11f]));
        maxy(hnd,1)=max(abs([bdens11; bdens11f]));
        %maxy=max(t1sd(1)*bdens11);miny=max(t1fsd(1)*bdens11f);
        %ylim([-1*miny maxy].*1.3)
        %xlim(th(1)+limspc.*max(abs([data1(:,1); data1f(:,1)]))*[-1 1])
        plot([min(abs(c11)),maxx(hnd,1)],[0,0],'k')
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
        ftopt=0;
        if ftopt==1
            ft=title(sprintf('\\color{blue}nested model fixed for \\nu=%0.1f',thf(2)),...
                'Interpreter','tex','FontWeight','normal','FontName','Helvetica',...
                'FontSize',afs);
        elseif ftopt==2
            if ~isinf(thf(2))
                ft=title(sprintf('nested model fixed for $$\\nu=%0.1f$$',thf(2)),...
                    'Interpreter','latex','FontSize',afs);
            else
                ft=title('nested model fixed for $$\\nu\rightarrow\infty$$',...
                    'Interpreter','latex','FontSize',afs);
            end
        end
        hold off

        % The smoothness estimates
        axes(ha((hnd-1)*numrows+3))
        % The infinity case will need a few changes
        hold on        
        kdelim=[min(data1(:,2)) max(data1(:,2))]+[-1 1].*datspc;    
        [a12,bdens12,c12]=kde(data1(:,2),2^kdep,kdelim(1),kdelim(2));         
        t12kp=bar(c12,bdens12,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        %kdelim=[min(data1f(:,2)) max(data1f(:,2))]+[-1 1].*1.1;       
        %[a12f,bdens12f,c12f]=kde(data1f(:,2),2^kdep,kdelim(1),kdelim(2)); 
        th02=scatter(th(2),0,ms,'k|'); 
        thf2=scatter(thhatf(2),0,floor(ms/4),'ko');
        if ~isinf(thf(2))
            thfix2=scatter(thf(2),0,floor(ms/2),'kd');
        end
        maxy(hnd,2)=max(abs([bdens12]));
        maxx(hnd,2)=max(abs([c12]));
        %maxy=max(t1sd(2)*bdens12);miny=0.1;                                 
        %ylim([-1*miny maxy].*1.3)                                            
        %if ~isinf(th(2))
        %    xlim(th(2)+limspc.*max(abs([data1(:,2); data1f(:,2)]))*[-1 1])  
        %end
        plot([min(abs(c12)),maxx(hnd,2)],[0,0],'k')
        if hnd==1                                                        
            ylabel('$\hat{\nu}$','Interpreter', 'latex')                 
            set(get(gca,'ylabel'),'rotation',0)                          
        end                                                                
        longticks
        set(gca,'YTickLabel',[])                                          
        %%% Here or later?
        if hnd~=3 % skip this for the nu->Inf case
            xtickvals1=gca().XAxis.TickValues;                                    
            nutickvals=[1/3 1/2 1 3/2 5/2];
            xtickvals=nutickvals;
            if any(xtickvals==th(2));
                nuticklabs={'$\frac{1}{3}$','$\frac{1}{2}$','$1$','$\frac{3}{2}$','$\frac{5}{2}$'};
                % for certain values, this can be difficult to read
                %nuticklabs{xtickvals==th(2)}=sprintf('$\\nu_0$=%0.2f',th(2));
            elseif ~isinf(th(2))
                [xtickvals,srtnu]=sort([xtickvals th(2)]);
                nuticklabs={'$\frac{1}{3}$','$\frac{1}{2}$','$1$','$\frac{3}{2}$','$\frac{5}{2}$'};
                % for certain values, this can be difficult to read
                % nuticklabs{srtnu==max(srtnu)}=sprintf('$\\nu_0$=%0.2f',th(2));
            else
                nuticklabs={'$\frac{1}{3}$','$\frac{1}{2}$','$1$','$\frac{3}{2}$','$\frac{5}{    2}$'};      
            end
            set(gca,'XTick',nutickvals,'XTickLabel',nuticklabs,'TickLabelInterpreter','latex','FontSize',afs);
        end
        box on                  
        hold off      

        % The range estimates
        axes(ha((hnd-1)*numrows+4))
        hold on                                                          
        kdelim=[min(data1(:,3)) max(data1(:,3))]+[-1 1].*datspc;          
        [a13,bdens13,c13]=kde(data1(:,3),2^kdep,kdelim(1),kdelim(2));        
        t13kp=bar(c13,bdens13,1,'FaceColor',clrs{1},'EdgeColor',clrs{1});
        kdelim=[min(data1f(:,3)) max(data1f(:,3))]+[-1 1].*datspc;    
        [a13f,bdens13f,c13f]=kde(data1f(:,3),2^kdep,kdelim(1),kdelim(2));
        t13fkp=bar(c13f,-1*bdens13f,1,'FaceColor',clrs{2},'EdgeColor',clrs{2});
        th03=scatter(th(3),0,ms,'k|');  
        thf3=scatter(thhatf(3),0,floor(ms/4),'ko');
        %maxy=max(t1sd(3)*bdens13);miny=max(t1fsd(3)*bdens13f);      
        %ylim([-1*miny maxy].*1.3)    
        %xlim(th(3)+limspc.*max(abs([data1(:,3); data1f(:,3)]))*[-1 1]) 
        maxy(hnd,3)=max(abs([bdens13; bdens13f]));
        maxx(hnd,3)=max(abs([c13 c13f]));
        plot([min(abs(c13)),maxx(hnd,3)],[0,0],'k')
        if hnd==1                                           
            ylabel('$\hat{\rho}$','Interpreter', 'latex')      
            set(get(gca,'ylabel'),'rotation',0)                            
        end                                     
        set(gca,'YTickLabel',[])               
        % xtickvals=gca().XAxis.TickValues;           
        % xtickvals=[xtickvals(1) xtickvals(end)];                          
        % if any(xtickvals==th(3));                                           
        %     rhticklabs=cellstr(num2str(xtickvals(:)./L,'%0.1g'));          
        %     rhticklabs{xtickvals==th(3)}=sprintf('$\\rho_0/[L]$=%0.1g',th(3)/L);                        
        % else                                                                
        %     [xtickvals,srtrh]=sort([xtickvals th(3)]);                      
        %     rhticklabs=cellstr(num2str(xtickvals(:)./L,'%0.1g'));           
        %     rhticklabs{srtrh==max(srtrh)}=sprintf('$\\rho_0/[L]$=%0.1g',th(3)/L);                       
        % end                                                              
        % set(gca,'XTick',xtickvals,'XTickLabel',rhticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');
        longticks
        box on                                                               
        hold off        
       %  %text(-maxdst/2,6,sprintf('$\\sigma^2_0=%0.1f$',th(1)),'Interpreter', 'latex');
       %  %text(-maxdst/2,4,sprintf('$\\nu_0=%0.1f$',th(2)),'Interpreter', 'latex');
       %  %text(-maxdst/2,2,sprintf('$\\rho_0=%0.1f$',th(3)),'Interpreter', 'latex');
       %  hold off

        axes(ha((hnd-1)*numrows+5))           
        hold on
        [h,cax,t]=plotit(Hx2,p2,[],th,0);
        longticks
        set(gca,'FontSize',9);
        if hnd>1;set(gca,'YTickLabels',[]);end
        box on;
        if ~isinf(thhatf(2))
           % title(sprintf('realization for $$\\hat{\\theta}=$$[%0.2f %0.2f %0.2f]',thhatf),...
           thhatf
           st2(hnd)=title(sprintf('$$\\hat{\\theta}=$$[%0.2f %0.2f %0.2f]',thhatf),...
                'Interpreter','latex','FontSize',afs);
        else
           % title(sprintf(...
           %     'realization for $$\\hat{\\theta}=$$[%0.2f $$\\infty$$ %0.2f]',thhatf(1),thhatf(3)),...
            st2(hnd)=title(sprintf(...
                '$$\\hat{\\theta}=$$[%0.2f $$\\infty$$ %0.2f]',thhatf(1),thhatf(3)),...
                'Interpreter','latex','FontSize',afs);
        end
        hold off

        % axes(ha((hnd-1)*numrows+1));rf1pos=gca().Position;
        % rf1npos=shrink(ha((hnd-1)*numrows+1),0.6,0.6);
        % % For some reason this isn't always effective, so the below lines are
        % % needed
        % if hnd>1
        %     ha((hnd-1)*numrows+1).Position(3)=ha(1).Position(3);
        %     ha((hnd-1)*numrows+1).Position(4)=ha(1).Position(4);
        %     ha((hnd-1)*numrows+1).Position(2)=ha(1).Position(2);
        %     %rf2npos=ha((hnd-1)*numrows+1).Position;
        % else
        %     ha((hnd-1)*numrows+1).Position(2)=rf1pos(2)-rf1npos(4)+rf1pos(4);
        % end

        % axes(ha((hnd-1)*numrows+5));rf2pos=gca().Position;
        % ha((hnd-1)*numrows+5).Position(3)=ha((hnd-1)*numrows+1).Position(3);
        % ha((hnd-1)*numrows+5).Position(4)=ha((hnd-1)*numrows+1).Position(4);
        % ha((hnd-1)*numrows+5).Position(1)=ha((hnd-1)*numrows+1).Position(1);
        % rf2npos=ha((hnd-1)*numrows+5).Position;
        % box on
        % %ha((hnd-1)*numrows+5).Position(2)=rf2pos(2)-rf2npos(4)+rf2pos(4);
        % 
        % thshr=2;
        % shrink(ha((hnd-1)*numrows+2),1,thshr)
        % shrink(ha((hnd-1)*numrows+3),1,thshr)
        % shrink(ha((hnd-1)*numrows+4),1,thshr)
        % movev(ha((hnd-1)*numrows+4),0.07)
        % %serre([ha((hnd-1)*numrows+3) ha((hnd-1)*numrows+4)],1/6,'down')
        % serre([ha((hnd-1)*numrows+2) ha((hnd-1)*numrows+3)],1/2,'down')
    end  
    % Now assign a common axis range for all of the kde subplots
    ylspc=1.2;
    stdOut=5;
    % Shared limits for finite smoothness panels
    ths=[th1; th2; th3];
    for hnd=1:3
        set(H(hnd+1,1),'XLim',ths(1,hnd)+max(t1sd(:,hnd)).*[-stdOut stdOut],...
            'YLim',max(maxy(1,hnd)).*[-ylspc ylspc]);
        set(H(hnd+1,2),'XLim',ths(2,hnd)+max(t1sd(:,hnd)).*[-stdOut stdOut],...
            'YLim',max(maxy(2,hnd)).*[-ylspc ylspc]);
    end
    % Treat labels for sigma and rho
    for hnd=1:3
    % select the sigma plots
        axes(H(2,hnd));
        xtickvals=gca().XAxis.TickValues;
        % For the infinity case, sigma2 values can grow large. Get rid of the 
        % first label so that the text has space
        if any(xtickvals==th(1));
            xticklabs=cellstr(num2str(xtickvals(:),'%0.1g'));
            xticklabs{xtickvals==th(1)}=sprintf('$\\sigma^2_0$=%0.1g',th(1));
            if max(xtickvals)~=th(1)
                xticklabs={xticklabs{1}; xticklabs{xtickvals==th(1)};...
                   xticklabs{end}};
                xtickvals=[xtickvals(1) th(1) xtickvals(end)];
            end
        else
            [xtickvals,srtsig]=sort([xtickvals th(1)]);
            xticklabs=cellstr(num2str(xtickvals(:),'%0.1g'));
            xticklabs{srtsig==max(srtsig)}=sprintf('$\\sigma^2_0$=%0.1g',th(1));
            if max(xtickvals)~=th(1)
               xticklabs={xticklabs{1}; xticklabs{srtsig==max(srtsig)};...
                    xticklabs{end}};
               xtickvals=[xtickvals(1) th(1) xtickvals(end)];
            end
        end
        set(gca,'XTick',xtickvals,'XTickLabel',xticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');

    % select the rho plots
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
        % For the infinity case, rho values can grow large. Get rid of the first
        % label so that the text has space
        if hnd==3
            xticklabs={0;xticklabs{2}; xticklabs{3}};
            xtickvals=[0 xtickvals(2) xtickvals(3)];
        end
        set(gca,'XTick',xtickvals,'XTickLabel',xticklabs,'TickLabelInterpreter','latex','FontSize',afs,'XTickLabelRotationMode','manual');
    end

    % Limits for the infinite smoothness panel
    keyboard
    for jnd=1:3
        %set(H(jnd+1,3),'XLim',[0 median(thhat_3(:,jnd))+std(thhat_3(:,jnd)).*stdOut],...
        set(H(jnd+1,3),'XLim',[0 maxx(3,jnd).*ylspc],...
            'YLim',max(maxy(3,jnd)).*[-1*ylspc ylspc]);
    end
    % Special treatment for nu x-axis for non-infinite case
    xlspc=1.1;cth=[th(1) 1.5 th(3)];
    set(H(3,1:2),'XLim',cth(2)+(max(abs(cth(2)-nutickvals))).*[-1*xlspc xlspc],'FontSize',afs,'XTickLabelRotationMode','manual');

    axes(H(1,1));rf1pos=gca().Position;axes(H(1,2));rf2pos=gca().Position;
    axes(H(1,3));rf3pos=gca().Position;
    axes(H(5,1));rf4pos=gca().Position;axes(H(5,2));rf5pos=gca().Position;
    axes(H(5,3));rf6pos=gca().Position;

    fshr=0.5;
    rf1npos=shrink(H(1,1),fshr,fshr);
    rf2npos=shrink(H(1,2),fshr,fshr);
    rf3npos=shrink(H(1,3),fshr,fshr);
    % For some reason SHRINK isn't always effective, so the below lines are
    % needed to have the latter fields have an equivalent size to the first
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

    thshr=1.7;
    shrink([H(2:4,:)],1,thshr)
    serre([H(2:4,1)],5/12,'down');serre([H(2:4,2)],5/12,'down');serre([H(2:4,3)],5/12,'down')
    movev(H(1,:),-0.02);
    movev(H,0.06);
    movev(H(5,:),-0.005);
    movev(H(2:4,:),-0.002);

    % Minimize space between the subplots
    serre([H(2:4,:)],1/3,'across')
    serre([H(2:4,:)],-3/64,'across')
    serre([H(1,:)'],-9/64,'across')
    serre([H(5,:)'],-9/64,'across')
    %H(1,1).XTickLabelRotation=0;H(1,2).XTickLabelRotation=0;H(1,3).XTickLabelRotation=0;
    axes(H(1,1));box on;axes(H(1,2));box on;axes(H(1,3));box on;
    axes(H(5,1));box on;axes(H(5,2));box on;axes(H(5,3));box on;
    for hnd=1:3
        st(hnd).Position(2)=st(hnd).Position(2)+2;
        st2(hnd).Position(2)=st2(hnd).Position(2)+2;
    end
    keyboard
% need to resize subfigures
    figdisp([],sprintf('specialfig5_%3.3i_%3.3i_%3.3i_fld_fix%i',p.NyNx,xver,fixflag),[],1);
end

function varargout=plotit(v,p,cr,th,mt,ah)
% Should use th(1) instead of halverange, shouldn't I...
vm=mean(v);
v=v-vm;
sox=std(v);
sax=[-2 2]*sox;
[h,cax]=imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],sax); axis xy image
%[h,cax]=imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],halverange(v,80,NaN)); axis ij image
%hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
hold on
hc=circ(pi*th(3),[],[p.NyNx(2)/5 p.NyNx(1)/5]);
delete(hc(2))
hc(1).Color='b';
xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
% Colorbar styling copied over from EGGERS2
ho2=getpos(ah);
[cb,xcb]=addcb([ho2(1)+ho2(3)*1.05  ho2(2)*1.07 ho2(3)/12 ho2(4)*0.75],...
  sax,sax,'kelicol',sox);
hold on
cbty=0.5+min(get(cb,'ylim'));
pgk=plot(get(cb,'xlim'),[cbty cbty],'k');
hold off
axes(cb)
xxb=2*indeks(get(cb,'xlim'),2);
yyb=get(cb,'ytick');
sl(1)=text(xxb,yyb(1),sprintf('%s','-2s'));
sl(2)=text(xxb,yyb(2),sprintf('%s','-1s'));
sl(3)=text(xxb,yyb(4),sprintf('%s','+1s'));
sl(4)=text(xxb,yyb(5),sprintf('%s','+2s'));
sl(5)=text(xxb,yyb(3),sprintf('%s','m'));
nolabels(cb,2)
delete(xcb)
defval('mt',1)
if mt==1
if ~isempty(th)
    if isstr(th)
        t=title(th);
    else
        t=title(sprintf('%i x %i | %s = [%g %g %g]',p.NyNx(1),p.NyNx(2),'\theta',...
                        th./[1 1 sqrt(prod(p.dydx))]));
    end
    movev(t,-p.NyNx(1)/25)
    varns={h,cax,t,cb};
end
else
    varns={h,cax,[],cb};
end
varargout=varns;
