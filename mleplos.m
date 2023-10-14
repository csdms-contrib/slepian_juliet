function varargout=mleplos(thhats,th0,covF0,covavhs,covXpix,E,v,params,name,thpix)
% MLEPLOS(thhats,th0,covF0,covavhs,covXpix,E,v,params,name,thpix)
%
% Graphical statistical evaluation of the maximum-likelihood inversion
% results from the suite MLEOS, MLEROS, MLEROS0, MLEOSL. Displays
% probability density functions of the estimates and makes
% quantile-quantile plots to ascertain normality. 
%
% INPUT:
%
% thhats     The estimated model parameter vector
% th0        The true model parameter vector
%            th0(1)=D    Isotropic flexural rigidity 
%            th0(2)=f2   The sub-surface to surface initial loading ratio 
%            th0(3)=r    The sub-surface to surface initial correlation coefficient
%            th0(3/4)=s2   The first Matern parameter, aka sigma^2 
%            th0(4/5)=nu   The second Matern parameter 
%            th0(5/6)=rho  The third Matern parameter 
% covF0      A (poor) covariance estimate based on the Fisher matrix at
%            the truth, which does not involve any of the data
% covavhs    The covariance matrix based on the median numerical Hessian
%            matrix near the individual estimates, from the diag file
% covXpix    A covariance matrix based on the numerical Hessian at a random estimate
% E          Young's modulus (not used for single fields)
% v          Poisson's ratio (not used for single fields)
% params     The structure with the fixed parameters from the experiment
% name       A name string for the title
% thpix      The example estimate, randomly picked up
%
% OUTPUT:
%
% ah,ha,yl,xl,tl  Various axis handles of the plots made
%
% EXAMPLE:
%
% This only gets used in MLEOS/MLEROS/MLEROS0/MLEOSL thus far, their 'demo2'
%
% Last modified by fjsimons-at-alum.mit.edu, 10/14/2023

defval('xver',1)

% Number of times the standard deviation for scale truncation
nstats=[-3:3]; fax=3;
sclth0=10.^round(log10(abs(th0)));
movit=0.01;
yls=[-0.0 0.75];
% Determines the rounding on the y axis 
rondo=1e2;
% Sets the format for the estimated/true plot labels

% The number of parameters
np=size(thhats,2);
if np==6
  labs={'D','f^2','r','\sigma^2','\nu','\rho'};
  labs0={'D_0','f^2_0','r_0','\sigma^2_0','\nu_0','\rho_0',};
  unts={'Nm' [] [] [] [] []};
elseif np==5
  labs={'D','f^2','\sigma^2','\nu','\rho',};
  labs0={'D_0','f^2_0','\sigma^2_0','\nu_0','\rho_0'};
  unts={'Nm' [] [] [] []};
elseif np==3
  labs={'\sigma^2','\nu','\rho',};
  labs0={'\sigma^2_0','\nu_0','\rho_0'};
  flabs={'variance \sigma^2','smoothness \nu','range \rho',};
  unts={[] [] []};
end

figure(1)
clf
[ah,ha]=krijetem(subnum(2,np));

disp(sprintf('\n'))

% For each of the parameters
for ind=1:np
  % The empirical means and standard deviations of the estimates
  mobs=nanmean(thhats(:,ind));
  sobs=nanstd(thhats(:,ind));
  % Collect them all
  mobss(ind)=nanmean(thhats(:,ind));
  sobss(ind)=nanstd(thhats(:,ind));

  % The means and standard deviations for any one estimate
  th0i=th0(ind);
  if ~isempty(covF0)
    % Error estimate based on the Fisher matrix at the truth
    stdF0=real(sqrt(covF0(ind,ind)));
  else
    stdF0=NaN;
  end
  % Collect them all
  stdF0s(ind)=stdF0;

  if ~isempty(covavhs)
    % Error estimate based on the median numerical Hessian matrix near estimate
    stdavhs=real(sqrt(covavhs(ind,ind)));
  else
    stdavhs=NaN;
  end
  % Collect them all
  stdavhss(ind)=stdavhs;

  if ~isempty(covXpix)
    % Error estimate based on one particular randomly picked numerical Hessian
    stdXpix=real(sqrt(covXpix(ind,ind)));
  else
      stdXpix=NaN;
  end
  % Collect them all
  stdXpixs(ind)=stdXpix;
  
  % HISTOGRAMS
  axes(ah(ind))
  % The "kernel density estimate"
  % Second input was a different default which we lowered
  [a,bdens,c]=kde(thhats(:,ind),2^8);
  if isinf(a) || any(bdens<-1e-10) || size(thhats,1)<50
    % If it's funky use the good old histogram
    [bdens,c]=hist(thhats(:,ind),max(size(thhats,1)/3,3));
    bdens=bdens/indeks(diff(c),1)/size(thhats(:,ind),1);
  end
  % This number is close to one... it's a proper density!
  if xver==1
    disp(sprintf('%s pdf normalization check by summation %g',...
                 upper(mfilename),sum(bdens)*indeks(diff(c),1)))
    disp(' ')
  end
  
  % Now plot it using a scale factor to remove the units from the y axis
  thhis(ind)=bar(c,sobs*bdens,1);
  set(ah(ind),'ylim',yls)
  % The markings
  stats=mobs+nstats*sobs;
  % What is the percentage captured within the range?
  nrx=20; nry=15;
  set(ah(ind),'XLim',stats([1 end]),'XTick',stats,'XTickLabel',nstats)

  % Truth and range based on stdF0
  hold on
  p0(ind)=plot([th0i th0i],ylim,'k-');
  halfup=indeks(ylim,1)+range(ylim)/2;
  ps(ind)=plot(th0i+[-1 1]*stdF0,...
	       [halfup halfup],'k-'); 
  % Didn't like this range bar in the end
  delete(ps(ind))
  
  % Estimate x-axis from observed means and variances
  xnorm=linspace(nstats(1),nstats(end),100)*sobs+mobs;
  % Normal distribution based on stdF0
  psF0(ind)=plot(xnorm,sobs*normpdf(xnorm,th0i,stdF0));
  % Based on the median numerical Hessian matrix at the estimates
  psavhs(ind)=plot(xnorm,sobs*normpdf(xnorm,th0i,stdavhs));
  % Based on one of them picked at random, numerical Hessian at estimate
  psXpix(ind)=plot(xnorm,sobs*normpdf(xnorm,th0i,stdXpix));
  % Based on the actually observed covariance of these data
  % In previous versions had used th0i/mobs here also, that didn't summarize it well
  pobs(ind)=plot(xnorm,sobs*normpdf(xnorm,th0i,sobs));

  % Some annotations
  % Experiment size, he giveth, then taketh away
  tu(ind)=text(stats(end)-range(stats)/nrx,indeks(ylim,2)-range(ylim)/nry,...
	      sprintf('N = %i',size(thhats(~isnan(sum(thhats,2))),1))); set(tu(ind),'horizon','r')
  fbb=fillbox(ext2lrtb(tu(ind),[],0.8),'w'); delete(tu(ind)); set(fbb,'EdgeColor','w')
  tu(ind)=text(stats(end)-range(stats)/nrx,indeks(ylim,2)-range(ylim)/nry,...
	      sprintf('N = %i',size(thhats(~isnan(sum(thhats,2))),1))); set(tu(ind),'horizon','r')

  % The percentage covered in the histogram that is being shown
  tt(ind)=text(stats(1)+range(stats)/nrx,indeks(ylim,2)-2*range(ylim)/nry,...
	      sprintf('s/%s = %5.2f','\sigma',sobs/stdavhs)); 
  fb=fillbox(ext2lrtb(tt(ind),[],0.8),'w'); delete(tt(ind)); set(fb,'EdgeColor','w')
  tt(ind)=text(stats(1)+range(stats)/nrx,indeks(ylim,2)-2*range(ylim)/nry,...
	      sprintf('s/%s = %5.2f','\sigma',sobs/stdavhs));

  % The ratio of the observed to the theoretical standard deviation
  t(ind)=text(stats(1)+range(stats)/nrx,indeks(ylim,2)-range(ylim)/nry,...
	      sprintf('%4.2f%s',...
		      sum(bdens([c>=stats(1)&c<=stats(end)])/sum(bdens)*100),...
		      '%'));
  hold off
  xl(ind)=xlabel(labs{ind});

  % QUANTILE-QUANTILE PLOTS  
  axes(ah(ind+np))
  h=qqplot(thhats(:,ind)); delete(h(2))
  set(h(1),'MarkerE','k')  
  set(h(3),'LineS','-','Color',grey)
  top(h(3),ah(ind+np))
  set(ah(ind+np),'xlim',nstats([1 end]),...
		'box','on','xtick',nstats,'XTickLabel',nstats)
  delete(get(ah(ind+np),'Ylabel'));
  delete(get(ah(ind+np),'Title'));
  delete(get(ah(ind+np),'XLabel'));

  set(ah(ind+np),'YLim',stats([1 end]),'YTickLabel',stats,...		       
		'YTickLabel',round(rondo*stats/sclth0(ind))/rondo);
  hold on
  e(ind)=plot(xlim,[mobs mobs],'k:');
  f(ind)=plot([0 0],ylim,'k:');
  bottom(e(ind),ah(ind+np))
  bottom(f(ind),ah(ind+np))
  set(ah(ind+np),'plotbox',[1 1 1])
  if sclth0(ind)~=1
    tl(ind)=title(sprintf('%s = %5.3f %s %4.0e %s',labs{ind},...
			  mobs/sclth0(ind),'x',...
			  sclth0(ind),unts{ind}));
    xl0(ind)=xlabel(sprintf('%s = %5.3f %s %4.0e %s',labs0{ind},...
			    th0(ind)/sclth0(ind),'x',...
			    sclth0(ind),unts{ind}));
  else
    tl(ind)=title(sprintf('%s = %5.3f %s',labs{ind},...
			  mobs/sclth0(ind),...
			  unts{ind}));
    xl0(ind)=xlabel(sprintf('%s = %5.3f %s',labs0{ind},...
			    th0(ind)/sclth0(ind),...
			    unts{ind}));
  end
  movev(xl0(ind),-range(ylim)/15)
  drawnow
end     

% Cosmetics
set(thhis(1:np),'FaceColor',grey,'EdgeColor',grey)
if np==6
  mv=0.125; mh=0.01; aps1=[0.8 1]; aps2=[1 1];
  set(thhis(4:6),'FaceColor',grey(9),'EdgeColor',grey(9))
elseif np==5
  mv=0.125; mh=0.01; aps1=[0.8 1]; aps2=[1 1];
  set(thhis(3:5),'FaceColor',grey(9),'EdgeColor',grey(9))
elseif np==3
  mv=0.1; mh=-0.075; aps1=[1.3 1]; aps2=[1.4 1.4];
end
for ind=1:np-1
  moveh(ha([1:2]+2*(ind-1)),(ind-np)*mh)
end
shrink(ah(1:np),aps1(1),aps1(2))
shrink(ah(np+1:end),aps2(1),aps2(2))

movev(ah(length(ah)/2+1:end),mv)
axes(ah(1))
yl=ylabel('posterior probability density');
axes(ah(4))
yl(2)=ylabel('sample values');
longticks(ah)
% Normal distribution based on stdF0
set(psF0,'linew',0.5,'color','k','LineS','--')
% Based on the median numerical Hessian matrix
set(psavhs,'linew',1.5,'color','k')
% Based on the randomly picked Hessian matrix
set(psXpix,'linew',0.5,'color','k')
% Based on the actually observed covariance of these data
set(pobs,'linew',1.5,'color',grey(3.5))

% Delete the one you know barely works
delete(psF0)

% Do this so the reduction looks slightly better
set(yl,'FontSize',12)
nolabels(ah(2:np),2)
%disp(sprintf('\n'))
%fig2print(gcf,'landscape')

% Stick the params here somewhere so we can continue to judge
movev(ah,-.1)
% If params isn't a structure, we're not in the right mindset
if isstruct(params)
  t=ostitle(ah,params,name); movev(t,.4)
end

% Here is the TRUTH and the COVF0 standard deviation
try 
    [answ,answs]=osansw(th0,covF0,E,v);
    disp(sprintf('%s',...
                 'Truth and Fisher-based covariance standard deviation'))
    disp(sprintf(answs,answ{:}))
    
    % Here is the RANDOMLY PICKED estimate and its NUMERICAL-HESSIAN based standard deviation
    [answ,answs]=osansw(thpix,covXpix,E,v);
    disp(sprintf('\n%s',...
                 'Example estimate and numerical-Hessian covariance standard deviation'))
    disp(sprintf(answs,answ{:}))

    % Here is the MEAN ESTIMATE and its OBSERVED-COVARIANCE-based standard
    % deviation - exactly like mobss and sobss==diag(sqrt(cov(thhats)))
    [answ,answs]=osansw(mean(thhats),cov(thhats),E,v);
    disp(sprintf('\n%s',...
                 'Mean estimate and ensemble-covariance standard deviation'))
    disp(sprintf(answs,answ{:}))

    % By the way, use THAT as a subtitle
    tt=supertit(ah(np+1:2*np),sprintf(answs,answ{:}));
end

if np>3; movev(tt,-4); else; movev(tt,-3.5); end

% Make basic x-y plots of the parameters
% SHOULD CALL THIS MLETHPLOS
if xver==1
    figure(2)
  clf
  pcomb=nchoosek(1:np,2);
  pstats=[-2 2]; tstats=[-3 3]; vstats=[-2 0 2];
  [ah,ha]=krijetem(subnum(1,3));
  
  % Scale everything
  mobss=mobss./sclth0;
  stdavhss=stdavhss./sclth0;
  sobss=sobss./sclth0;
  thhats=thhats./repmat(sclth0,size(thhats,1),1);
  th0=th0./sclth0;
  
  for ind=1:np
    axes(ah(ind))
    % Find the pairwise combinations
    p1=pcomb(ind,1); p2=pcomb(ind,2);

    % Observed means and theoretical standard deviations
    t1(ind)=plot(mobss(p1)+pstats*stdavhss(p1),...
		 [mobss(p2) mobss(p2)]); hold on
    t2(ind)=plot([mobss(p1) mobss(p1)],...
		 mobss(p2)+pstats*stdavhss(p2));
    set([t1 t2],'Color',grey)
    % The parameter estimates
    p(ind)=plot(thhats(:,p1),thhats(:,p2),'o'); 

    % Observed means and observed standard deviations
    m(ind)=plot(mobss(p1),mobss(p2),'v');
    o1(ind)=plot(mobss(p1)+pstats*sobss(p1),...
		[mobss(p2) mobss(p2)],'LineWidth',2);
    o2(ind)=plot([mobss(p1) mobss(p1)],...
		mobss(p2)+pstats*sobss(p2),'LineWidth',2);
    hold off
    % Truths
    try
      set(ah(ind),'xtick',round(100*[th0(p1)+vstats*sobss(p1)])/100,...
		  'ytick',round(100*[th0(p2)+vstats*sobss(p2)])/100)
    end
    axis square;  grid on
    xlim(th0(p1)+tstats*sobss(p1))
    ylim(th0(p2)+tstats*sobss(p2))
    % Color mix
    cmix=[0 0 0]; cmix([p1 p2])=1/2;
    set([p(ind) m(ind)],'MarkerFaceColor',cmix,'MarkerEdgeColor',cmix,'MarkerSize',2)
    % Cosmetix
    delete([o1(ind) o2(ind)])
    xlabel(flabs{p1})
    ylabel(flabs{p2})
  end
  longticks(ah)
  %seemax([ah(1) ah(2)],1)
  %seemax([ah(2) ah(3)],2)
  titi=ostitle(ah,params,name); movev(titi,-2)
  try
      tt=supertit(ah(1:np),sprintf(answs,answ{:})); movev(tt,-7)
  end
end

% Output
varns={ah,ha,yl,xl,tl};
varargout=varns(1:nargout);

% Subfunction to compute standard deviations
