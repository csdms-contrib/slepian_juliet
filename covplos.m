function varargout=covplos(oneortwo,sclcovX,obscov,params,thhats,E,v,cblox,ifinv)
% COVPLOS(oneortwo,sclcovX,obscov,params,thhats,E,v,cblox,ifinv)
%
% Makes a nice covariance plot of the maximum-likelihood estimation results
%
% INPUT:
%
% oneortwo  1 Plots only prediction
%           2 Plots prediction and observation [default]
% sclcovX   The scaled covariance estimate
% obscov    The scaled sample covariance matrix
% params    The parameters from the experiment, for the title only
% thhats    The estimated model parameter vector, whose length we need
% E         Young's modulus (not used for single fields and only for title)
% v         Poisson's ratio (not used for single fields and only for title)
% cblox     Is the colorbar 'hor' or 'ver' [default]
% ifinv      Ordered inversion flags for [s2 nu rho], e.g. [1 0 1]; if [1 1 1],
%            the loglikelihood contours will be estimated and included with the
%            second figure, otherwise, the appropriate subplots will be removed
%
% OUTPUT:
%
% ah,yl     Axis handles to what was just plotted
%
% Last modified by fjsimons-at-alum.mit.edu, 05/20/2025
% Last modified by olwalbert-at-princeton.edu, 05/20/2025

defval('oneortwo',2)
defval('cblox','ver')
defval('ifinv',[1 1 1])

% Number of paramters
np=length(sclcovX);

% Calculate proper covariance of scores, how about using the average estimate
% maybe should bring that in line with MLEPLOS which uses the truth
thhat=mean(thhats,1);
try
    % DFMTX
    cvg=covgammiosl(thhat,params,2);
catch
    % gradient sample
    cvg=covgammiosl(thhat,params,1);
end
% Calculate true covariance of estimated parameters and replace what came in
sclcovX=covthosl(thhat,params,cvg,[1 1 1]);
sclcovX=sclcovX./[diag(sqrt(sclcovX))*diag(sqrt(sclcovX))'];

% Print this out
obscov
sclcovX

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
  im=imagesc(obscov);
  im.AlphaData=~isnan(im.CData);
  covannotate(ah(2),np,ifinv)
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
movev(t,.5);
set(t,'FontSize',12);
moveh(t,-0.1);
% We are quoting the TRUTHS and the 1xstandard deviation based on COVX
% Here is the MEAN ESTIMATE and its OBSERVED-COVARIANCE-based standard deviation
[answ,answs]=osansw(mean(thhats),cov(thhats),E,v);
disp(sprintf('\n%s',...
             'Mean estimate and ensemble-covariance standard deviation'))
disp(sprintf(answs,answ{:}))
% By the way, use THAT as a title
tt=supertit(ah(1:2),sprintf(answs,answ{:}));
movev(tt,-4.5)
set(tt,'FontSize',12);
moveh(tt,-0.1);

% Output
varns={ah,yl};
varargout=varns(1:nargout);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function covannotate(ah,np,ifinv)

defval('ifinv',repelem(1,np))

if np==6
  xc=[3.5 3.5];
  yc=[0.5 6.5];
  labs={'D','f^2','r','\sigma^2','\nu','\rho'};
elseif np==5
  xc=[2.5 2.5];
  yc=[0.5 5.5];
  labs={'D','f^2','r','\sigma^2','\nu','\rho'};
elseif np==3
  labs={'\sigma^2','\nu','\rho'};
end
labs(~ifinv)={''};

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
set(ah,'xtick',1:np,'XTickLabel',labs,'FontSize',12)
set(ah,'ytick',1:np,'YTickLabel',labs,'FontSize',12)
longticks(ah)
movev(ah,.05)

