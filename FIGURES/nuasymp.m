function nuasymp(th,p,dth,xver)
% NUASYMP(TH,P,DTH,XVER)
%
% Creates a plot showing the isotropic, stationary Matern spectrum as a radial
% projection and its partial derivative with respect to a Matern parameter
% selected by the dth index in context of the wave-vector and the parameter
%
% Input:
%
% th    Isotropic Matern spectral parameters
% p     Grid parameters p.NyNx and p.dydx for calculating
%       the wavevector, and p.blurs for calculating the spectrum
% dth   Matern parameter index of the parameter that will be allowed to vary, 
%       1 :: s2
%       2 :: nu
%       3 :: rh
% xver  how many curves to plot for each function? default to 3
%
% Output:
%
% <>    a figure will be saved locally
%
% Examples:
%
% p=[];p.NyNx=[32 32];p.dydx=[1 1];nuasymp([1 1 5],p,1,3)
%
% Last modified by olwalbert@princeton.edu 08/30/2024

defval('th',[1e6 1.5 2e4])
defval('xver',3)
defval('dth',2)
defval('p',[])

% Set axes to semilogx or loglog? Flag 1 or 0
ll=1;

% How many dimensions? 
dim=2;

if isempty(p)
    p.NyNx=[128 128];
    p.dydx=[1e3 1e3];
    p.blurs=Inf;
end

% For wavevector dependent subfigure, define a dense wavevector and devise
% indices for xver values of th(dth) that are of interest; for parameter 
% dependent subfigure, define a dense parameter vector for th(dth) and devise 
% indices for wavenumbers of interest 

% Wavenumber axes, 2D axes flattened to 1D radial, inspired from EGGERS2
mfromkm=1000;
Nx=p.NyNx(2);
dx=p.dydx(2);
% Wavenumber limits [rad/m]
[k,~,~,dci]=knum2([1 Nx],[1 (Nx-1)*dx]);
k=k(dci(2):end);
% Finer grid for S(k)
k=unique([linspace(k(1),k(2),50)...
    linspace(k(2),k(3),50)...
    linspace(k(3),k(4),50)...
    k(4:end)]);
% Spatial axes, for distance [m]
x=linspace(0,(Nx-1)*dx,100);
% New range and ticks for Ky
newr=[-0.1 1.1];
newy=[0 1/3 2/3 1];
% Major tick wavelengths [km]
xtil=[100 1000 10000];
% Wavenumber limits [rad/km]
klims=[k(2) max(k)]*mfromkm;

% Labels in order of dth indices
labs={'\sigma^2','\nu','\rho'};
axlabs={'variance, \sigma^2','smoothness, \nu','range, \rho'};

% Matern parameters of interest 
thalld={[th(1)/2:th(1)/50:th(1)*3],[0.05:0.1:20],[th(3)/2:th(3)/50:th(3)*3]};
thsub={[2*th(1)/3 th(1) 2*th(1)],[1/2 3/2 5/2],[2*th(3)/3 th(3) th(3)*2]};
thall{dth}=[thalld{dth} thsub{dth}];
% The fluctuation scale [m^4], following from EGGERS2; will need to be
% recalculated if dth==1 or dth==3
switch dth
    case 1
        S0d=thalld{dth}.*pi.*th(3).^2/4;
        S0=thall{dth}.*pi.*th(3).^2/4;
    case 2
        S0d=th(1).*pi.*th(3).^2/4;
        S0=S0d;
    case 3
        S0d=th(1).*pi.*thalld{dth}.^2/4;
        S0=th(1).*pi.*thall{dth}.^2/4;
end
[thall{dth},temp]=sort(thall{dth});[~,thind]=sort(temp);
thind=thind(end-xver+1:end);
N0=length(thalld{dth});
N1=length(thsub{dth});
N2=length(thall{dth});

% Wavenumbers of interest 
knd=[floor(length(k)/20); floor(length(k)/9); floor(length(k)/5)];
ksub=k(knd)';

% Preallocate space for the periodograms, their derivatives, and their finite
% difference approximations 
Sks=zeros(N0,length(k));
Sk1s=zeros(N1,length(k));
Sk2s=zeros(N0,length(ksub));
dSk1s=zeros(N1,length(k));
dSk2s=zeros(N0,length(ksub));
cdiffSks=zeros(length(k),N2);
% Preallocate space for the L2 norms of ratios between periodograms of 
% increasing th
rSk=zeros(N2-1,1);
rdSk=zeros(N2-1,1);

% Calculate expected periodogram for the dense th vector for all wavenumbers
mnd=1;
for ind=1:N2
    switch dth
        case 1
            th1=[thall{dth}(ind) th(2) th(3)];
        case 2
            th1=[th(1) thall{dth}(ind) th(3)];
        case 3
            th1=[th(1) th(2) thall{dth}(ind)];
    end
    Sks(ind,:)=maternos(k,th1,[],dim,1); 
    if any(ind==thind)
        % if the current Matern parameter is one that we will plot profiles for,
        % take the opportunity to calculate the partial derivative and store the
        % periodogram separately
        dSk1s(mnd,:)=maternos(k,th1,dth,dim,1);
        Sk1s(mnd,:)=Sks(ind,:);
        mnd=mnd+1;
    else
        %Sks(ind-mnd+1,:)=maternos(k,th1,[],dim,1);
        dSk2s(ind-mnd+1,:)=maternos(ksub,th1,dth,dim,1); 
        if ind>1
            r=Sks(ind,:)./Sks(ind-1,:)-1;
            %r=r(~isnan(r));r=r(~isinf(r));
            rSk(ind-1)=norm(r,2);
            rd=dSk2s(ind-mnd+1,:)./dSk2s(ind-mnd,:)-1; 
            %rd=rd(~isnan(rd));rd=rd(~isinf(rd));
            rdSk(ind-1)=norm(rd,2);
        end
    end
end

% Let's normalize here
Sks=(Sks'./S0)';
switch dth
    case 1 
        Sk1s=(Sk1s'./S0(thind))';
        dSk1s=(dSk1s'./S0(thind))';
    case 2
        Sk1s=(Sk1s'/S0)';
        dSk1s=(dSk1s'/S0)';
    case 3
        Sk1s=(Sk1s'./S0(thind))';
        dSk1s=(dSk1s'./S0(thind))';
end
Sk2s=(Sk2s'./S0d)';
dSk2s=(dSk2s'./S0d)';

% Parse dense Sks for wavenumbers of interest for the parameter dependent plots
% and calculate the partial derivative
for jnd=1:xver
    temp=Sks(:,knd(jnd));
    Sk2s(:,jnd)=[temp(1:thind(1)-1); temp(thind(1)+1:thind(2)-1);...
                temp(thind(2)+1:thind(3)-1); temp(thind(3)+1:end)];
                %Sks(:,knd(jnd)); %maternos(ksub,th2,[],dim,1);
end
% Since we skip the thind indices when calculating the periodogram norm ratios,
% lets set any zeros in the ratios to NaNs for plotting
rSk(rSk==0)=NaN;
rSk(find(isnan(rSk)==1)+1)=NaN;
rdSk(rdSk==0)=NaN;

% Calculate the central difference approximation for wavenumber dependence
for ind=1:length(k)
    for jnd=2:N2-1
        % across columns, jnd: th(dth); down rows, ind: k
        % want central difference for changes in th(dth)
        cdiffSks(ind,jnd)=(Sks(jnd+1,ind)-Sks(jnd-1,ind))/(2*(thall{dth}(jnd+1)-thall{dth}(jnd-1)));
    end
% Forward difference for first obs 
    cdiffSks(ind,1)=(Sks(2,ind)-Sks(1,ind))/(thall{dth}(2)-thall{dth}(1));
% Backward difference for last obs 
    cdiffSks(ind,N2)=(Sks(N2,ind)-Sks(N2-1,ind))/(thall{dth}(N2)-thall{dth}(N2-1));
end
% We only need the central difference approximation for columns indexed
% by thind and rows indexed by knd
cdiffSk1s=cdiffSks(:,thind);
cdiffSk2s=cdiffSks(knd,:);

% The colors, marker sizes
clrs=colormap(parula(4));
ms=4;

% Create the figure with three column plots
figure(1);clf;
[ah,ha,H]=krijetem(subnum(1,3));

% Create the first suplot
axes(ha(1))
hold on
for ind=1:xver
    plot([ksub(ind) ksub(ind)],[-0.5 1.5],...
         ':', 'Color',[0.65 0.65 0.65],'LineWidth',0.1);
end
if ll
    pl1=semilogx(k,Sk1s,'-',...
                 k,dSk1s,'-.',...
                 k,2*cdiffSk1s,'o');
    pl1(end-2).MarkerSize=ms;
    pl1(end-1).MarkerSize=ms;
    pl1(end  ).MarkerSize=ms;
else
    plot(k,Sk1s,'-',k,dSk1s,':',k,cdiffSk1s,'o');
end
xl(1)=xlabel('k, wavenumber (rad/km)');
yl(1)=ylabel('normalized by S0');

% Let's create a concise legend where the colors (values of parameter) and
% marker style (function) only appear once each
for ind=1:xver
   % Specify the xver different values by color
   LC1(ind)=plot(nan,nan,...
            'LineWidth',2,'DisplayName',...
            sprintf('$%s_%i$=%0.2f',labs{dth},ind,thsub{dth}(ind)));
end
% Specify the functions by line or marker style
LC1(xver+1)=plot(nan,nan,...
            'LineWidth',2,'Color','k','DisplayName',...
            sprintf('$S(k; %s_i)$',labs{dth}));
LC1(xver+2)=plot(nan,nan,...
            'LineStyle','-.','LineWidth',2,'Color','k','DisplayName',...
            sprintf('$\\partial S (k;%s_i) / \\partial %s$',...
                    labs{dth},labs{dth}));
LC1(xver+3)=plot(nan,nan,...
            'ko','DisplayName',...
            sprintf('$\\tilde{\\partial S}(k;%s_i) / \\partial %s$',...
                    labs{dth},labs{dth}));
LCgnd1=legend(LC1,'Location','best','NumColumns',1,'Interpreter','latex');
LCgnd1.ItemTokenSize=[12 12];

% Clean up the axes
yscl=0.025;
xlim([k(1) max(k)]);
ylim([min(min(dSk1s))-yscl max(max(Sk1s))+yscl])
if ll
    xscale('log');
end
box on;
set(gca,"ColorOrder",clrs(1:3,:))
ytck=get(gca,'YTick');set(gca,'YTick',ytck(1:2:end));
longticks;
hold off

% Create the middle subplot
axes(ha(2))
hold on
for ind=1:xver
    vl2(ind)=plot([thsub{dth}(ind) thsub{dth}(ind)],[1e-10 1e5],...
             ':','Color',[0.65 0.65 0.65],'LineWidth',0.5);
end
if ll
    pl2=loglog(thalld{dth},Sk2s,'-',...
               thalld{dth},dSk2s,'-.',...
               thall{dth}(2:end-1),2*cdiffSk2s(:,2:end-1),'o');
    pl2(end-2).MarkerSize=ms;
    pl2(end-1).MarkerSize=ms;
    pl2(end  ).MarkerSize=ms;
else
    plot(thall{dth},Sk2s,'-',...
         thall{dth},dSk2s,'-.',...
         thall{dth}(2:end-1),2*cdiffSk2s(2:end-1,:),'o');
end

% Let's create a concise legend where the colors (values of parameter) and
% marker style (function) only appear once each
for ind=1:xver
   % Specify the xver different values by color
   LC2(ind)=plot(nan,nan,...
            'LineWidth',2,'DisplayName',...
            sprintf('$k_%i$=%0.1e',ind,ksub(ind)));
end
% Specify the functions by marker/line style
LC2(xver+1)=plot(nan,nan,...
            'LineWidth',2,'Color','k','DisplayName',...
            sprintf('$S(%s;k_i)$',labs{dth}));
LC2(xver+2)=plot(nan,nan,...
            'LineWidth',2,'LineStyle','-.','Color','k','DisplayName',...
            sprintf('$\\partial S (%s;k_i) / \\partial %s$',...
            labs{dth},labs{dth}));
LC2(xver+3)=plot(nan,nan,...
            'ko','DisplayName',...
            sprintf('$\\tilde{\\partial S}(%s;k_i) / \\partial %s$',...
            labs{dth},labs{dth}));
LCgnd2=legend(LC2,'Location','best','NumColumns',1,'Interpreter','latex');
LCgnd2.ItemTokenSize=[12 12];

% Clean up the axes
xlabel(axlabs{dth})
yscl=[2 1.1];
ylim([min(min([dSk2s Sk2s])) max(max([dSk2s Sk2s]))].*yscl)
if ll
    xscale('log');
    yscale('log');
end
box on;
set(gca,"ColorOrder",clrs(1:3,:))
ytck=get(gca,'YTick');set(gca,'YTick',ytck(1:2:end));
longticks;
hold off

% Make a title with information relevant to creating the periodograms
switch dth
    case 1
        tistr=sprintf('$\\nu=%0.1g,\\rho=%0.1g$',th(2),th(3));
    case 2
        tistr=sprintf('$\\sigma^2 = %0.1g, \\rho = %0.1g$',th(1),th(3));
    case 3
        tistr=sprintf('$\\sigma^2 = %0.1g, \\nu = %0.1g$',th(1),th(2));
end
ti=title(sprintf('%s;  %ix%i grid;  %ix%i km',...
         tistr,p.NyNx,p.NyNx.*p.dydx/1e3),'Interpreter','latex');
movev(ti,0.01);

% Make the final subplot
axes(ha(3))
hold on
for ind=1:xver
    vl3(ind)=plot([thsub{dth}(ind) thsub{dth}(ind)],[1e-20 max(max([rSk; rdSk]))+1],...
             ':','Color',[0.65 0.65 0.65],'LineWidth',0.1);
end
set(gca,"ColorOrder",...
    [0 0.4470 0.7410;0.9290 0.6940 0.1250;0.6350 0.0780 0.180])
if ll
    pl3=loglog(thall{dth}(2:end),rSk,'-',...
                 thall{dth}(2:end),rdSk,'-.');
else
    plot(thall{dth}(2:end),rSk,'-',...
         thall{dth}(2:end),rdSk,'-.');
end

% Let's create a nice legend 
LC3(1)=plot(nan,nan,...
            'LineWidth',2,'DisplayName',...
            sprintf('$||\\frac{S(k; %s_{i+1})}{S(k; %s_{i})} - 1||_2$',...
                    labs{dth},labs{dth}));
LC3(2)=plot(nan,nan,...
            '-.','LineWidth',2,'DisplayName',...
   sprintf('$||\\frac{ \\partial S(k; %s_{i+1})/ \\partial %s}{ \\partial S(k;%s{i} / \\partial %s } - 1 ||_2$',...
             labs{dth},labs{dth},labs{dth},labs{dth}));
LCgnd3=legend(LC3,'NumColumns',1,'Interpreter','latex');
LCgnd3.ItemTokenSize=[12 12];

% Clean up the axes
xlabel(axlabs{dth})
yscl=[0.9 1.1];
ylim([min(min([rSk rdSk])) max(max([rSk rdSk]))].*yscl)
if ll
    xscale('log');
    yscale('log');
end
box on;
set(gca,"ColorOrder",clrs(1,:))
longticks;
hold off

% Save the figure
switch dth
    case 1
        tith=[th(2) th(3)];
    case 2
        tith=[th(1) th(3)];
    case 3
        tith=[th(1) th(2)];
end
figdisp([],sprintf('nuasymp_%3.3i_%3.3i_%idth_%i_%i',p.NyNx,dth,tith),[],1);
