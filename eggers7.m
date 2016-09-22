function eggers7
% EGGERS7
%
% Makes FIGURE 4 of Olhede et al. (2017), illustrating the empirical
% performance of estimating the parameters of an isotropic Matern process,
% as a function of its parameters (variance, differentiability, range) and
% for a suite of field sizes. The data plotted in this figure are generated
% by invoking SIMULOSL('demo4').
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 09/19/2016

varos={'\sigma','\nu','\rho'};

% Metric conversion
mfromkm=1000;

% Matern parameters, in standard units
% The parameters for the top row in the published work
th0t=[1e6 2.0 5e4];
% The parameters for the bottom row in the published work
th0b=[1e6 3.0 10e4];

% Alternative set... used for EGGERS4, save as -final-alt
%th0t=[1e6 2.5 2e4];
%th0b=[1e6 2.5 4e4];

% Experimental parameters
fields={'dydx','NyNx','blurs','kiso','quart'};
defstruct('params',fields,{[10 10]*1e3,[64 64],-1,NaN,0})
% Maximum field size and discretization of the sizes tried
N=128; rindj=[2:2:N];

% Number of simulations per processor, higher is better for small grids
npr=5;
% Number of processors, must agree with your machine and SIMULOSL
NumWorkers=8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the figure: four panels with different parameters
clf; [ah,ha,H]=krijetem(subnum(2,3));

% All maximum-likelihood cases
for pix=1:3
  axes(ah(pix))
  [~,~,ax(pix)]=...
      simulosl('demo4',th0t,params,'mle',N,rindj,npr,ah,[],[],pix);
  delete(ax(pix))
  % Automatic gain wasn't too good looking, so try again
  ylim([-th0t(pix)/10 th0t(pix)*2.5])
  set(ah(pix),'ytick',[0:0.5:2]*th0t(pix),'yticklabel',{'0' '' '1' '' '2'})
end
% Meaningful title
tl(1)=supertit(ah(2),sprintf('%s = %3.1f km ; %s = %3.1f ; %s = %i km',...
			'\sigma',sqrt(th0t(1))/mfromkm,...
			'\nu',th0t(2),...
			'\rho',th0t(3)/mfromkm));

% All maximum-likelihood cases
for pix=1:3
  axes(ah(pix+3))
  [~,~,ax(pix+3)]=...
      simulosl('demo4',th0b,params,'mle',N,rindj,npr,ah,[],[],pix);
  delete(ax(pix+3))
  ylim([-th0b(pix)/10 th0b(pix)*2.5])
  set(ah(pix+3),'ytick',[0:0.5:2]*th0b(pix),'yticklabel',{'0' '' '1' '' '2'})
end
% Meaningful title
tl(2)=supertit(ah(5),sprintf('%s = %3.1f km ; %s = %3.1f ; %s = %i km',...
			'\sigma',sqrt(th0b(1))/mfromkm,...
			'\nu',th0b(2),...
			'\rho',th0b(3)/mfromkm));

% Cosmetics to match EGGERS1 and EGGERS4 as best we can
for ind=1:length(ah); yl(ind)=get(ah(ind),'ylabel'); end
for ind=1:length(ah); xl(ind)=get(ah(ind),'xlabel'); end

% If it goes full-width in the LaTeX
set(tl,'FontSize',8)
set([ah xl yl],'FontS',7)
set(findobj('Marker','o'),'Marker','.','LineStyle','none')

delete(findobj('Color','blue'))
delete(findobj('Color','magenta'))

set(ah,'ygrid','on','box','on')
set(yl(1),'string','estimates relative to truth')
set(yl(4),'string','estimates relative to truth')
set(xl([4 5 6]),'string',sprintf('grid size (%s%s)','\pi','\rho'))
nolabels(ha([3 4 5 6]),2); delete(yl([2 3 5 6]))
nolabels(ah([1 2 3]),1); delete(xl([1 2 3]))
shrink(ah,1,2)
serre(H,[],'across'); serre(H',0.7,'down')
pitx=0:2:8;
pibx=0:4;
set(ah(1:3),'xtick',pi*pitx*th0t(3)/mfromkm,'xticklabel',pitx)
set(ah(4:6),'xtick',pi*pibx*th0b(3)/mfromkm,'xticklabel',pibx)

for pix=1:3
  axes(ah(pix))
  [bh(pix),th(pix)]=boxtex('ur',ah(pix),varos{pix},8,1,1,1.3,1);
end

for pix=1:3
  axes(ah(pix+3))
  [bh(pix+3),th(pix+3)]=boxtex('ur',ah(pix+3),varos{pix},8,1,1,1.3,1);
end

% Delete the black vertical markers
set(intersect(intersect(findobj('LineWidth',0.5),findobj('Color','k')),findobj('MarkerSize',6)),...
    'LineStyle',':')
% Set the data-size-limited boundary on the correlation length, see SIMULOSL
dl=findobj('Color','r'); set(dl,'LineWidth',0.25,'Color','k')
bottom(dl,ah(3))
bottom(dl,ah(6))

% Set the medians in there also, see SIMULOSL
dm=findobj('Color','c'); set(dm,'Marker','+','MarkerSize',2,'Color','k','LineStyle','none')
delete(dm)

% Set the trimmed-variances in there also, see SIMULOSL
dv=findobj('Color','y'); set(dv,'Marker','none','Color','k','LineStyle','-')
delete(dv)

movev(tl(1),-2.9); 
movev(tl(2),-0.8)
moveh(tl,0.15)

% Print to file
fig2print(gcf,'portrait')
figdisp([],[],[],2,'epsc','epstopdf');
