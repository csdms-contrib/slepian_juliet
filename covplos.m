function varargout=covplos(oneortwo,sclcovX,obscov,covX,params,thhats,th0,E,v,cblox)
% COVPLOS(oneortwo,sclcovX,obscov,covX,params,thhats,th0,E,v,cblox)
%
% Makes a nice covariance plot of the MLE results
%
% INPUT:
%
% oneortwo  1 Plots only prediction
%           2 Plots prediction and observation [default]
% sclcovX   The scaled predicted covariance matrix
% obscov    The scaled observed covariance matrix
% covX      The predicted covariance matrix before scaling (used for the title)
% params    The parameters from the experiment, for the title only
% thhats    The estimated model parameter vector, whose length we need
% th0       The true model parameter vector
%           th0(1)=D    Isotropic flexural rigidity 
%           th0(2)=f2   The sub-surface to surface initial loading ratio 
%           th0(3)=r    The sub-surface to surface initial correlation coefficient
%           th0(1/4)=s2   The first Matern parameter, aka sigma^2 
%           th0(2/5)=nu   The second Matern parameter 
%           th0(3/6)=rho  The third Matern parameter 
% E         Young's modulus (not used for single fields and only for title)
% v         Poisson's ratio (not used for single fields and only for title)
% cblox     Is the colorbar 'hor' or 'ver' [default]
%
% OUTPUT:
%
% ah,yl     Axis handles to what was just plotted
%
% Last modified by fjsimons-at-alum.mit.edu, 10/17/2017

defval('oneortwo',2)
defval('cblox','ver')

% Number of paramters
np=length(sclcovX);

switch oneortwo
 case 1
  % Plot
  clf
  fig2print(gcf,'portrait')  
  ah=gca;
  imagesc(sclcovX)
  covannotate(ah,np)

  shrink(ah,1.5,1.5)
  cb=colorbar('hor');
  shrink(cb,1/1.225,1.5)
  movev(cb,-.175); longticks(cb,2)
  axes(cb)
  xlabel('normalized predicted covariance matrix')
  set(cb,'ydir','rev')
 case 2
  clf
  fig2print(gcf,'portrait')  
  ah=krijetem(subnum(2,2));
  delete(ah(3:4))
  ah=ah(1:2);

  axes(ah(1))
  imagesc(sclcovX)
  covannotate(ah(1),np)
  t(1)=title('predicted');

  axes(ah(2))
  imagesc(obscov)
  covannotate(ah(2),np)
  t(2)=title('observed');

  set(ah(2),'yaxisl','r')
  movev(ah,-.25)
  serre(ah,1,'across')
  
  switch cblox
   case 'hor'
    cb=colorbarf('hor',get(gca,'FontSize'),get(gca,'FontName'),...
		 [0.3 0.3 0.3 0.02]);
    axes(cb)
    xlabel('normalized variance/covariance matrix')
    moveh(cb,0.125)
    movev(cb,-0.02)
   case 'ver'
    cb=colorbarf('ver',get(gca,'FontSize'),get(gca,'FontName'),...
		 [0.175 0.385 0.02 getpos(ah(1),4)]);
    axes(cb)
    yl=ylabel('normalized covariance matrix');
    set(cb,'yaxisloc','l')
  end
  longticks(cb,2)
end

% Stick the params here somewhere so we can continue to judge
t=ostitle(ah,params,[],length(thhats(:,1)));
movev(t,.5)
% We are quoting the TRUTHS and the 1xstandard deviation based on COVX
[answ,answs]=osansw(th0,covX,E,v);

tt=supertit(ah,sprintf(answs,answ{:}));
movev(tt,-5)

moveh(tt,-0.075)

% Output
varns={ah,yl};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function covannotate(ah,np)

if np==6
  xc=[3.5 3.5];
  yc=[0.5 6.5];
  labs={'D','f2','r','s2','nu','rho'};
elseif np==5
  xc=[2.5 2.5];
  yc=[0.5 5.5];
  labs={'D','f2','r','s2','nu','rho'};
elseif np==3
  labs={'s2','nu','rho'};
end

kelicol
caxis([-1 1])

if np>3
  hold on
  p(1)=plot(xc,yc,'k');
  p(2)=plot(yc,xc,'k');
  set(p,'Linew',2)
  hold off
end
axis image

% Labels and such
set(ah,'xtick',1:np,'xtickl',labs)
set(ah,'ytick',1:np,'ytickl',labs)
longticks(ah)
movev(ah,.05)

