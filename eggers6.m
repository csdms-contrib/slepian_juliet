function eggers6
% EGGERS6
%
% Makes FIGURE 5 of Olhede et al. (2017), illustrating the chi-squaredness
% of the residuals and the normality of the test statistic for
% chi-squaredness. May have to run it twice for proper proportions.
%
% Last modified by fjsimons-at-alum.mit.edu, 07/07/2015

% Set the experimental parameters
th0=[1e6 2.5 2e4];
fields={'dydx','NyNx','blurs','kiso','quart'};
defstruct('params',fields,{[10 10]*1e3,[64 64],2,NaN,0});

% Do one SIMULOSL simulation and MLEOSL recovery
[Hx,th0,params,k1,Hk1,Sb,Lb]=simulosl(th0,params);
[thhat,~,~,scl,~,p,Hk2,k2]=mleosl(Hx,[],params,[],[],[]);

% Checks and balances... note that Hk2 has been demeaned
diferm(Hk1(~~k1),Hk2(~~k1),5)
diferm(k1,k2)

% Must make some decent axes
clf
ah=krijetem(subnum(2,3));

% Plot these guys bjutifooly
[cb,xl,t,tt]=mlechiplos(4,Hk1,thhat,scl,params,ah(1:3));

keyboard

delete(t)

% Then here plot the results from the series of tests in a new figure
% that we pillage and kill afterwards
figure
whatscreated='16-Jun-2015-64-2';
[cobs,mobs,nobs,th0,p,momx]=mleosl('demo2',whatscreated);
close(gcf)
% Now plot what we reaped, i.e. momx(:,3)

% HISTOGRAM PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(4))
mx3=momx(:,3);
% Take out the Inf which may occur at zero wavenumber
[bdens,c]=hist(mx3,4*round(log(length(mx3))));
bdens=bdens/indeks(diff(c),1)/length(mx3);
bb=bar(c,bdens,1);
varibal='s_X^2';
% Each of the below should be df/2
t(1)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
		   varibal,nanmean(mx3),...
		   varibal,nanvar(mx3)));
set(t,'FontS',9)
set(bb,'FaceC',grey)
hold on
sfax=4;
varpred=8/length(k1(~~k1));
xll=[1-sfax*sqrt(varpred) 1+sfax*sqrt(varpred)];
xls=linspace(1-sfax*sqrt(varpred),1+sfax*sqrt(varpred),sfax+1);
refs=linspace(xll(1),xll(2),100);
plot(refs,normpdf(refs,1,sqrt(varpred)),'Linew',1,'Color','k')
axis square
xlim(xll)
set(ah(4),'xtick',xls,'xtickl',round(xls*100)/100)
longticks(ah(4))
xl(1)=xlabel(varibal); 
yl(1)=ylabel('probability density');
ylim([0 11])
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save this position as it gets messed up afterwards
posh=getpos(ah(4));

% QUANTILE-QUANTILE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(5))
h=qqplot(mx3,ProbDistUnivParam('normal',[1 sqrt(varpred)]));
axis image; box on
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
set(ah(5),'xtick',xls,'xtickl',xlc,...
          'ytick',xls,'ytickl',xlc)
longticks(ah(5))
xl(2)=xlabel(sprintf('observed %s',varibal));
yl(2)=ylabel(sprintf('predicted %s',varibal));
axis square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The SPACE-DOMAIN DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(6))
lx=sqrt(length(Hx));
xlis=[0.5 lx/2 lx+0.5];
sfax=2;
cax=[-sfax sfax]*sqrt(th0(1))/1000;
imagefnan([xlis(1) xlis(end)],[xlis(end) xlis(1)],...
	  reshape(Hx-mean(Hx),params.NyNx)/1000,...
	  'gray',cax,[],1); axis image
xstr='simulated field [km]';
cborien='vert';
set(ah(6),'xtick',xlis,'xtickl',[-lx/2 0 lx/2],...
	  'ytick',xlis,'ytickl',[-lx/2 0 lx/2])


[cb(2),xcb]=addcb(cborien,cax,cax,'gray',sqrt(th0(1))/1000,1);
set(xcb,'string',xstr)
set([xcb cb(2)],'fonts',12)
set(ah(6),'position',...
 	  [getpos(ah(6),1) getpos(ah(6),2) getpos(ah(5),3) ...
 	   getpos(ah(5),4)])

axes(ah(6))
xl(3)=xlabel('position index'); 
yl(3)=ylabel('position index'); 
moveh(yl(3),params.NyNx(1)/30)
longticks(ah(6))
axes(cb(2))
moveh(cb(2),.075)
set(cb(2),'yaxisl','r')
shrink(cb(2),1.3,1)
set(cb(2),'position',...
       [getpos(cb(2),1) getpos(ah(6),2) getpos(cb(2),3) getpos(ah(6),4)])
shrink(cb(2),0.95,1.27)
moveh(xcb,15)
moveh([ah(4:6) cb(2)],-.025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last minute cosmetics
for ind=4:6
  set(ah(ind),'position',getpos(ah(ind)).*[1 1 0 0]+[0 0 posh(3) posh(4)]*0.99)
end
moveh(ah(4:6),[getpos(ah(1),1)-posh(1)*0.9])

movev([ah(4:6) cb(2)],-0.01)
for ind=1:3
  axes(ah(3+ind))
  yls=ylim; xls=xlim;
  movev(xl(ind),-range(xls)/25)
  if ind==2
    moveh(yl(ind),range(xls)/25)
  end
end
movev(xl(3),+range(yls)/25/2)
moveh(cb(2),.005)
moveh(xcb,-5)

% Write out
figna=figdisp([],[],[],1);
system(sprintf('epstopdf %s.eps',figna));

