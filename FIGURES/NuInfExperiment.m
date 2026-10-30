function NuInfExperiment(numvals,numreals,algo)
% For several realizations with increasing values of nu, how well can we recover
% nu? At what point can we effectively treat nu as infinite?
%INPUT
%
% numvals     how many (logarithmically spaced) values of nu between 0 and 100 
% numreals    how many realizations
% algo        optimization algorithm, defaults to unconstrained
%
% EXAMPLE
%
% NuInfExperiment(20,48,'unc')
%
% Last modified by olwalbert-at-princeton.edu, 04/21/2025

defval('numvals',20)
defval('numreals',48)
defval('algo','unc')

savef=1;
fnam='NuInfExp_';
fnam=append(fnam,"_optalgo_",algo,'_numreals_',num2str(numreals));

p=[];p.dydx=[1 1];p.NyNx=[128 128];p.blurs=Inf;

s2=1; rh=3;
nus=logspace(0,2,numvals);

try
    % Do we already data saved?
    filename=append(fnam,'sthhat_est3params.mat');
    ests=importdata(filename);
    disp('loading previously saved data');
catch
    thhats=zeros(numvals,numreals,3);
    scls=zeros(numvals,numreals,3);
    
    % first case: estimate all three parameters
    for ind=1:numvals
        th0=[s2 nus(ind) rh];
        thhat=cell(1,numreals);
        scl=cell(1,numreals);
        Hx=cell(prod(p.NyNx),numreals);
        thhati=[];scli=[];
        while size(thhati,1)<numreals
            parfor jj=1:numreals
                Hx{jj}=simulosl(th0,p);
                [thhat{jj},~,~,scl{jj}]=mleosl(Hx{jj},[],p,algo,[],[],1);
            end
            thhati=[thhati; cat(1,thhat{:})];scli=[scli; cat(1,scl{:})];
            idxs=isnan(thhati(:,1));
            thhati(idxs,:)=[];scli(idxs,:)=[];
        end
        thhats(ind,:,:)=thhati(1:numreals,:);scls(ind,:,:)=scli(1:numreals,:); 
        clear Hx;
        %while isnan(thhat)
        %    try
        %        [thhat,~,~,scl]=mleosl(Hx,[],p,algo,[],[],1); 
        %    catch
        %        keyboard
        %    end
        %end
        %thhats(ind,:)=thhat;scls(ind,:)=scl;
    end
    ests=thhats.*scls;

    if savef==1
        save(append(fnam,'sthhat_est3params'),'ests');
    end
    clear thhat scl thhati scli;
end

try
    filename=append(fnam,'sthhat_est2params.mat');
    est2=importdata(filename);
    disp('loading previously saved data');
catch
    % second case: fix nu to progressively larger values and estimate 2 parameters
    thhats2=zeros(numvals,numreals,3);
    scls2=zeros(numvals,numreals,3);
    try 
        for ind=1:numvals
            th0=[s2 nus(ind) rh];
            thhat=cell(1,numreals);
            scl=cell(1,numreals);
            Hx=cell(prod(p.NyNx),numreals);
            thhati=[];scli=[];
            while size(thhati,1)<numreals   
                parfor jj=1:numreals
                    Hx{jj}=simulosl(th0,p);
                    [thhat{jj},~,~,scl{jj}]=mleosl(Hx{jj},[],p,algo,[],th0,2);
               end
               thhati=[thhati; cat(1,thhat{:})];scli=[scli; cat(1,scl{:})];
               idxs=isnan(thhati(:,1));
               thhati(idxs,:)=[];scli(idxs,:)=[];
            end
            thhats2(ind,:,:)=thhati(1:numreals,:);scls2(ind,:,:)=scli(1:numreals,:);
            clear Hx;
        end
    catch
        for ind=1:numvals
            th0=[s2 nus(ind) rh];
            thhat=cell(1,numreals);
            scl=cell(1,numreals);
            Hx=cell(prod(p.NyNx),numreals);
            thhati=[];scli=[];
            while size(thhati,1)<numreals   
                parfor jj=1:numreals
                    Hx{jj}=simulosl(th0,p);
                    [thhat{jj},~,~,scl{jj}]=mleosl(Hx{jj},[],p,algo,[],th0,2);
               end
               thhati=[thhati; cat(1,thhat{:})];scli=[scli; cat(1,scl{:})];
               idxs=isnan(thhati(:,1));
               thhati(idxs,:)=[];scli(idxs,:)=[];
            end
            thhats2(ind,:,:)=thhati(1:numreals,:);scls2(ind,:,:)=scli(1:numreals,:);
            clear Hx;
        end
    end
    est2=thhats2.*scls2;
    if savef==1
        save(append(fnam,'sthhat_est2params'),'est2');
    end
    clear thhat scl thhati scli;
end

try
    filename=append(fnam,'sthhat_est2paramsInf.mat');
    est3=importdata(filename);
    disp('loading previously saved data');
catch
    % third case: fix nu to infinity, estimate 2 parameters
    thhats3=zeros(numreals,3);
    scls3=zeros(numreals,3);
    
    th0=[s2 Inf rh];
    thhat=cell(1,numreals);
    scl=cell(1,numreals);
    Hx=cell(prod(p.NyNx),numreals);
    thhati=[];scli=[];
    while size(thhati,1)<numreals
        parfor jj=1:numreals
            Hx{jj}=simulosl(th0,p);
            [thhat{jj},~,~,scl{jj}]=mleosl(Hx{jj},[],p,algo,[],th0,2);
        end
        thhati=[thhati; cat(1,thhat{:})];scli=[scli; cat(1,scl{:})];
        idxs=isnan(thhati(:,1));
        thhati(idxs,:)=[];scli(idxs,:)=[];
    end
    thhats3=thhati(1:numreals,:);scls3=scli(1:numreals,:);
    est3=thhats3.*scls3;
    if savef==1
        save(append(fnam,'sthhat_est2paramsInf'),'est3');
    end
end

mhat1=permute(mean(ests,2),[1,3,2]);sehat1=permute(std(ests,[],2)/sqrt(numreals),[1,3,2]);
mhat2=permute(mean(est2,2),[1,3,2]);sehat2=permute(std(est2,[],2)/sqrt(numreals),[1,3,2]);
mhat3=mean(est3);sehat3=std(est3)/sqrt(numreals);

% color friendly palette
clrs=["#1B9E77","#D95F02","#7570B3"];
ms=16;
figure();                                                                   
% publication figure flg==1; original is not
flg=1;
if flg==1
  [ah,ha]=krijetem(subnum(1,2));
  axes(ah(1))
  hold on
  errorbar(nus,mhat1(:,1),sehat1(:,1),'^','MarkerEdgeColor',clrs{1},...
      'MarkerFaceColor',clrs{1},'Color',clrs{1},...
      'DisplayName','$$\hat{\sigma}^2$$')
  errorbar(nus,mhat1(:,2),sehat1(:,2),'s','MarkerEdgeColor',clrs{2},...
      'MarkerFaceColor',clrs{2},'Color',clrs{2},'DisplayName','$$\hat{\nu}$$')
  errorbar(nus,mhat1(:,3),sehat1(:,3),'o','MarkerEdgeColor',clrs{3},...
      'MarkerFaceColor',clrs{3},'Color',clrs{3},'DisplayName','$$\hat{\rho}$$')
  yline(s2,':','Color',clrs{1},'LineWidth',2,'DisplayName','$$\sigma^2_0$$')
  plot(nus,nus,'LineStyle',':','Color',clrs{2},'LineWidth',2,'DisplayName','$$\nu_0$$')
  yline(rh,':','Color',clrs{3},'LineWidth',2,'DisplayName','$$\rho_0$$')
  grid on
  set(gca,'xscale','log','yscale','log','ylim',[0.89 105],...
          'xlim',[0.89 115],'xtick',[1 10 100])
  box on
  longticks
  legend('Location','northwest','Interpreter','latex','NumColumns',1)
  xlabel('\nu_0');
  ylabel('$$\hat{\theta}=[\hat{\sigma}^2, \hat{\nu},\hat{\rho}]$$',...
      'Interpreter','latex');

  axes(ah(2))
  hold on
  errorbar(nus,mhat2(:,1),sehat2(:,1),'^','MarkerEdgeColor',clrs{1},...
      'MarkerFaceColor',clrs{1},'Color',clrs{1},'DisplayName','$$\hat{\sigma}^2$$')
  errorbar(nus,mhat2(:,2),sehat2(:,2),'s','MarkerEdgeColor',clrs{2},...
      'MarkerFaceColor',clrs{2},'Color',clrs{2},'DisplayName','$$\hat{\nu}$$')
  errorbar(nus,mhat2(:,3),sehat2(:,3),'o','MarkerEdgeColor',clrs{3},...
      'MarkerFaceColor',clrs{3},'Color',clrs{3},...
      'DisplayName','$$\hat{\rho}$$')
  yline(s2,':','Color',clrs{1},'LineWidth',2,'DisplayName','$$\sigma^2_0$$')
  plot(nus,nus,'LineStyle',':','Color',clrs{2},'LineWidth',2,'DisplayName','$$\nu_0$$')
  yline(rh,':','Color',clrs{3},'LineWidth',2,'DisplayName','$$\rho_0$$')
  grid on
  set(gca,'xscale','log','yscale','log','ylim',[0.89 105],...
          'xlim',[0.89 115],'xtick',[1 10 100])
  box on
  longticks
  xlabel('\nu_0');
  ylabel('$$\hat{\theta}=[\hat{\sigma}^2, \nu_0,\hat{\rho}]$$',...
      'Interpreter','latex');
  h=title(sprintf(...
      'Estimation of Matern parameters using "%s" optimizer on a %ix%i grid',...
      algo,p.NyNx),'Interpreter','tex');
  set(ah(1),'FontSize',12)
  set(ah(2),'FontSize',12)
  movev(h,0.05);
  moveh(ah(1),-0.01)
  moveh(ah(2),0.01)
else
  [ah,ha]=krijetem(subnum(1,3));                                              
  axes(ah(1))
  hold on
  %scatter(nus,ests(:,1),ms,'Marker','^','MarkerEdgeColor',clrs{1},...
  errorbar(nus,mhat1(:,1),sehat1(:,1),'^','MarkerEdgeColor',clrs{1},...
      'MarkerFaceColor',clrs{1},'Color',clrs{1},...
      'DisplayName','$$\hat{\sigma}^2$$')
  %scatter(nus,ests(:,2),ms,'Marker','s','MarkerEdgeColor',clrs{2},...
  errorbar(nus,mhat1(:,2),sehat1(:,2),'s','MarkerEdgeColor',clrs{2},...
      'MarkerFaceColor',clrs{2},'Color',clrs{2},'DisplayName','$$\hat{\nu}$$')
  %scatter(nus,ests(:,3),ms,'MarkerEdgeColor',clrs{3},...
  errorbar(nus,mhat1(:,3),sehat1(:,3),'o','MarkerEdgeColor',clrs{3},...
      'MarkerFaceColor',clrs{3},'Color',clrs{3},'DisplayName','$$\hat{\rho}$$')
  yline(s2,':','Color',clrs{1},'LineWidth',2,'DisplayName','$$\sigma^2_0$$')
  plot(nus,nus,'LineStyle',':','Color',clrs{2},'LineWidth',2,'DisplayName','$$\nu_0$$')
  yline(rh,':','Color',clrs{3},'LineWidth',2,'DisplayName','$$\rho_0$$')
  grid on
  set(gca,'xscale','log','yscale','log','ylim',[0.1 150])
  box on
  legend('Location','northwest','Interpreter','latex','NumColumns',1)
  xlabel('\nu_0');
  ylabel('$$\hat{\theta}=[\hat{\sigma}^2, \hat{\nu},\hat{\rho}]$$',...
      'Interpreter','latex');
  hold off
  
  axes(ah(2))
  hold on
  %scatter(nus,est2(:,1),ms,'Marker','^','MarkerEdgeColor',clrs{1},...
  errorbar(nus,mhat2(:,1),sehat2(:,1),'^','MarkerEdgeColor',clrs{1},...
      'MarkerFaceColor',clrs{1},'Color',clrs{1},'DisplayName','$$\hat{\sigma}^2$$')
  %scatter(nus,est2(:,2),ms,'Marker','s','MarkerEdgeColor',clrs{2},...
  errorbar(nus,mhat2(:,2),sehat2(:,2),'s','MarkerEdgeColor',clrs{2},...
      'MarkerFaceColor',clrs{2},'Color',clrs{2},'DisplayName','$$\hat{\nu}$$')
  %scatter(nus,est2(:,3),ms,'MarkerEdgeColor',clrs{3},'MarkerFaceColor',clrs{3},...
  errorbar(nus,mhat2(:,3),sehat2(:,3),'o','MarkerEdgeColor',clrs{3},...
      'MarkerFaceColor',clrs{3},'Color',clrs{3},...
      'DisplayName','$$\hat{\rho}$$')
  yline(s2,':','Color',clrs{1},'LineWidth',2,'DisplayName','$$\sigma^2_0$$')
  plot(nus,nus,'LineStyle',':','Color',clrs{2},'LineWidth',2,'DisplayName','$$\nu_0$$')
  yline(rh,':','Color',clrs{3},'LineWidth',2,'DisplayName','$$\rho_0$$')
  grid on
  set(gca,'xscale','log','yscale','log','ylim',[0.1 150])
  box on
  xlabel('\nu_0');
  ylabel('$$\hat{\theta}=[\hat{\sigma}^2, \nu_0,\hat{\rho}]$$',...
      'Interpreter','latex');
  h=title(sprintf(...
      'Estimation of Matern parameters using "%s" optimizer on a %ix%i grid',...
      algo,p.NyNx),'Interpreter','tex');
  movev(h,0.05);
  hold off
  
  axes(ah(3))
  hold on
  %scatter(ones(size(est3(:,1))),est3(:,1),ms,'Marker','^','MarkerEdgeColor',...
  errorbar(ones(size(est3(:,1))),mhat3(:,1),sehat3(:,1),'^',...
      'MarkerEdgeColor',clrs{1},'Color',clrs{1},...
      'MarkerFaceColor',clrs{1},'DisplayName','$$\hat{\sigma}^2$$')
  %scatter(ones(size(est3(:,3))),est3(:,3),ms,'MarkerEdgeColor',clrs{3},...
  errorbar(ones(size(est3(:,3))),mhat3(:,3),sehat3(:,3),'o',...
      'MarkerEdgeColor',clrs{3},'Color',clrs{3},...
      'MarkerFaceColor',clrs{3},'DisplayName','$$\hat{\rho}$$')
  yline(s2,':','Color',clrs{1},'LineWidth',2,'DisplayName','$$\sigma^2_0$$')
  yline(rh,':','Color',clrs{3},'LineWidth',2,'DisplayName','$$\rho_0$$')
  grid on
  set(gca,'xscale','linear','yscale','log','ylim',[0.1 150])
  box on
  xticklabels({'','\infty'});
  xlabel('\nu_0');
  ylabel('$$\hat{\theta}=[\hat{\sigma}^2, \nu\to\infty,\hat{\rho}]$$',...
      'Interpreter','latex');
  
  % Clean up spacing
  axes(ah(1))
  moveh(gca,-0.05);
  axes(ah(3))
  moveh(gca,0.05);
end
keyboard
saveas(gcf,append('NuInfExperiment_',string(date)),'epsc')
