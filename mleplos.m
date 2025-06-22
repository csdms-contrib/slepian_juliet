function varargout=mleplos(thhats,th0,covF0,covavhs,covXpix,E,v,params,name,thpix,ifinv)
% MLEPLOS(thhats,th0,covF0,covavhs,covXpix,E,v,params,name,thpix,ifinv)
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
% ifinv      Ordered inversion flags for [s2 nu rho], e.g. [1 0 1]; if [1 1 1],
%            the loglikelihood contours will be estimated and included with the
%            second figure, otherwise, the appropriate subplots will be removed
%
% OUTPUT:
%
% ah,ha,yl,xl,tl,p,m,o1,o2,ep,ec,axlm,cvth,ttt  Various axis handles of the plots made
%
% This gets used in 'demo2' of MLEOS/MLEROS/MLEROS0/MLEOSL and MASKIT/MUCKIT
%
% Tested on 8.5.0.197613 (R2015a)
%
% Last modified by olwalbert-at-princeton.edu, 05/30/2025
% Last modified by fjsimons-at-alum.mit.edu, 05/30/2025

defval('xver',1)
defval('ifinv',[1 1 1])

% Number of times the standard deviation for scale truncation
nstats=[-3:3]; fax=3;
pstats=[-1 1];
tstats=[-3 3];
vstats=[-2 -1 0 1 2];
sclth0=10.^round(log10(abs(th0)));
movit=0.01;
yls=[-0.0 0.6];
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

% Append the scaling of each Matern parameter axis if it is not 1
for i=1:np                                                       
    if sclth0(i)~=1
        try
            flabs{i}=append(flabs{i},sprintf(' (x 10^{%d})',log10(sclth0(i))));
        catch
            % Maybe this isn't necessary for R2015a
            flabs{i}=sprintf('%s %s',flabs{i},sprintf(' (x 10^{%d})',log10(sclth0(i))));
        end
    end
end  
% Set the taper
if isfield(params,'mask')
  if strcmp(params.mask,'random')
    mcl=str2double(name(1:2))*0.01;
    [~,I]=muckit(randn(params.NyNx),params,mcl);
    params.taper=I;
  else
    % Calculate the analytical covariance for the Matern parameter, grid, and taper
    [~,I]=maskit(randn(params.NyNx),params);
    if strcmp(name,params.mask)
      params.taper=I;
    else
      % Assume this is a call for the anti-mask from MASKIT('DEMO2')
      params.taper=~I;
    end
  end
else
  params.taper=1;
end

% Calculate covariance of scores
try
    % DFMTX -- may fail when the grid is too large
    cvg=covgammiosl(th0,params,2);
catch
    try 
        % gradient sample -- may fail if positive definite embedding cannot be
        % found
        cvg=covgammiosl(th0,params,1);
    catch
        % diagonals -- only the last resort because it is time consuming
        cvg=covgammiosl(th0,params,3);
    end
end

% Calculate true covariance of estimated parameters
% The mean estimate doesn't exist, especially given the tradeoffs
%cvth=covthosl(nanmean(thhats,1),params,cvg,[1 1 1]);
% Somewhat oddly variable, this appears
%cvth=covthosl(thhats(randi(size(thhats,1)),:),params,cvg,[1 1 1]);
cvth=covthosl(th0,params,cvg,[1 1 1]);

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
      % If you have nothing maybe make it zero so the second plot still goes
      stdavhs=0;
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

  % The standard deviation from the calculated covariance
  scth=sqrt(cvth(ind,ind));
  % Collect them all
  scths(ind)=scth;
  
  % HISTOGRAMS
  axes(ah(ind))
  % The "kernel density estimate"
  % Second input was a different default which we lowered
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
  p0(ind)=plot([th0i th0i],ylim,'k:');
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
  % Based on the actually observed mean and covariance of these data
  pobs(ind)=plot(xnorm,sobs*normpdf(xnorm,mobs,sobs));
  % Include the analytically calculated covariance 
  pclc(ind)=plot(xnorm,sobs*normpdf(xnorm,mobs,scth));

  % Some annotations
  % Experiment size, he giveth, then taketh away
  tu(ind)=text(stats(end)-range(stats)/nrx,indeks(ylim,2)-range(ylim)/nry,...
	      sprintf('N = %i',size(thhats(~isnan(sum(thhats,2))),1))); set(tu(ind),'horizon','r')
  fbb=fillbox(ext2lrtb(tu(ind),[],0.8),'w'); delete(tu(ind)); set(fbb,'EdgeColor','w')
  tu(ind)=text(stats(end)-range(stats)/nrx,indeks(ylim,2)-range(ylim)/nry,...
	      sprintf('N = %i',size(thhats(~isnan(sum(thhats,2))),1))); set(tu(ind),'horizon','r')

  % The ratio of the observed to the theoretical standard deviation
  tt(ind)=text(stats(1)+range(stats)/nrx,indeks(ylim,2)-2*range(ylim)/nry,...
	       sprintf('s/%s = %5.2f','\sigma',sobs/scth));
  fb=fillbox(ext2lrtb(tt(ind),1.20,1),'w'); delete(tt(ind)); set(fb,'EdgeColor','w')
  tt(ind)=text(stats(1)+range(stats)/nrx,indeks(ylim,2)-2*range(ylim)/nry,...
	      sprintf('s/%s = %5.2f','\sigma',sobs/scth));

  % The percentage covered in the histogram that is being shown
  t(ind)=text(stats(1)+range(stats)/nrx,indeks(ylim,2)-range(ylim)/nry,...
	      sprintf('%4.2f%s',...
		      sum(bdens([c>=stats(1)&c<=stats(end)])/sum(bdens)*100),...
		      '%'));
  hold off
  xl(ind)=xlabel(labs{ind});

  % QUANTILE-QUANTILE PLOTS  
  axes(ah(ind+np))
  h=qqplot(thhats(:,ind)); delete(h(2))
  set(h(1),'MarkerEdgeColor','k')  
  set(h(3),'LineStyle','-','Color',grey)
  top(h(3),ah(ind+np))
  set(ah(ind+np),'xlim',nstats([1 end]),...
		'box','on','Xtick',nstats,'XTickLabel',nstats)
  delete(get(ah(ind+np),'Ylabel'));
  delete(get(ah(ind+np),'Title'));
  delete(get(ah(ind+np),'XLabel'));
  % Label the sample values
  set(ah(ind+np),'YLim',stats([1 end]),'YTick',stats,...		       
     'YTickLabel',round(rondo*stats/sclth0(ind))/rondo);
  hold on
  % Plot the sample mean and two standard deviations
  %e(ind)=plot(xlim,[mobs mobs],'-','Color',grey);
  e{ind}=plot(xlim,repmat(stats([2 4 6]),2,1),'-','Color',grey);
  f{ind}=plot(repmat([-2 0 2],2,1),ylim,'k:');
  for jnd=1:length(e{ind})
      bottom(e{ind}(jnd),ah(ind+np))
      bottom(f{ind}(jnd),ah(ind+np))
  end
  set(ah(ind+np),'plotbox',[1 1 1])
  if sclth0(ind)~=1
    tl(ind)=title(sprintf('%s = %5.3f %s %4.0e %s',labs{ind},...
			  mobs/sclth0(ind),'x',...
			  sclth0(ind),unts{ind}));
    set(tl(ind),'FontWeight','normal');
    xl0(ind)=xlabel(sprintf('%s = %5.3f %s %4.0e %s',labs0{ind},...
			    th0(ind)/sclth0(ind),'x',...
			    sclth0(ind),unts{ind}));
  else
    tl(ind)=title(sprintf('%s = %5.3f %s',labs{ind},...
			  mobs/sclth0(ind),...
			  unts{ind}));
    set(tl(ind),'FontWeight','normal');
    xl0(ind)=xlabel(sprintf('%s = %5.3f %s',labs0{ind},...
			    th0(ind)/sclth0(ind),...
			    unts{ind}));
  end
  movev(tl(ind),range(ylim)/20)
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
  % mv=0.05; mh=-0.06; aps1=[1.3 1]; aps2=[1.3 1.3];
end
for ind=1:np-1
  moveh(ha([1:2]+2*(ind-1)),(ind-np)*mh)
end
shrink(ah(1:np),aps1(1),aps1(2))
shrink(ah(np+1:end),aps2(1),aps2(2))

movev(ah(length(ah)/2+1:end),mv)

axes(ah(1))
yl=ylabel('probability density');
axes(ah(4))
yl(2)=ylabel('sample values');
longticks(ah)
% Normal distribution based on stdF0
set(psF0,'LineWidth',0.5,'Color','k','LineStyle','--')
% Based on the median numerical Hessian matrix
set(psavhs,'LineWidth',1.5,'Color','k')
% Based on the randomly picked Hessian matrix
set(psXpix,'LineWidth',0.5,'Color','k')
% Based on the actually observed covariance of these data
set(pobs,'LineWidth',1.5,'Color',grey(3.5))
% The latest prediction
set(pclc,'LineWidth',1,'Color','k')
% Since 2025 turns out the Fisher matrix IS being blurred with -1
delete(psF0)
% Since 2025 we know the Hessians/Fisher are no good for blurred with Inf
delete(psavhs)
delete(psXpix)

% Do this so the reduction looks slightly better
% set(yl,'FontSize',12)
nolabels(ah(2:np),2)
%disp(sprintf('\n'))
%fig2print(gcf,'landscape')

% Stick the params here somewhere so we can continue to judge
movev(ah,-.1)
% If params isn't a structure, we're not in the right mindset
if isstruct(params)
  t=ostitle(ah,params,name,length(thhats(:,1))); movev(t,.7);
  set(t,'FontSize',12);
end

try 
    % Here is the TRUTH and the COVF0 standard deviation
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
    [answ,answs,pm]=osansw(mean(thhats),cov(thhats),E,v);
    disp(sprintf('\n%s',...
                 'Mean estimate and ensemble-covariance standard deviation'))
    disp(sprintf(answs,answ{:}))

    % By the way, use THAT as a subtitle (i.e., the last one in this list!)
    tt=supertit(ah(np+1:2*np),sprintf('%s\n%s%s',sprintf(answs,answ{:}),pm,...
                                      'one sigma uncertainty based on the ensemble'));
end

% Late-breaking
if np>3; movev(tt,-4); else; movev(tt,-3.5); end
set(tt,'FontSize',10)

if any(ifinv~=[1 1 1])
    delete(ah(~repmat(ifinv,1,2)))
end

% Make basic x-y plots of the parameters
% FJS SHOULD CALL THIS MLEXPLOS
if xver==1
    figure(2)
    fig2print(gcf,'landscape')
    clf
    pcomb=nchoosek(1:np,2);
    [ah,ha]=krijetem(subnum(1,3));

    % Scale everything
    mobss=mobss./sclth0;
    stdavhss=stdavhss./sclth0;
    sobss=sobss./sclth0;
    thhats=thhats./repmat(sclth0,size(thhats,1),1);
    th0=th0./sclth0;

    % Plot error ellipses without using ERROR_ELLIPSE
    % https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
    defval('cl',0.68)
    % Check this for three variables
    % s=-2*log(1-cl);
    s=chi2inv(cl,2);
    
    t=linspace(0,2*pi);

    for ind=1:np
        axes(ah(ind))
        % Find the pairwise combinations for the cross-plot convention: s2-nu,
        % s2-rho, nu-rho
        p1=pcomb(ind,1); p2=pcomb(ind,2);

        % Observed means and AVERAGE NUMERICAL HESSIAN standard deviations
        t1(ind)=plot(mobss(p1)+pstats*stdavhss(p1),...
		     [mobss(p2) mobss(p2)]); hold on
        t2(ind)=plot([mobss(p1) mobss(p1)],...
		     mobss(p2)+pstats*stdavhss(p2));
        set([t1(ind) t2(ind)],'Color',grey)
        % Observed means and EXACT standard deviations which should be great
        c1(ind)=plot(mobss(p1)+pstats*scths(p1),...
		     [mobss(p2) mobss(p2)]);
        c2(ind)=plot([mobss(p1) mobss(p1)],...
		     mobss(p2)+pstats*scths(p2));
        set([c1(ind) c2(ind)],'Color',grey)

        % The parameter estimates
        p(ind)=plot(thhats(:,p1),thhats(:,p2),'o'); 

        % Observed means and OBSERVED standard deviations
        m(ind)=plot(mobss(p1),mobss(p2),'v');
        o1(ind)=plot(mobss(p1)+pstats*sobss(p1),...
		     [mobss(p2) mobss(p2)],'LineWidth',1);
        o2(ind)=plot([mobss(p1) mobss(p1)],...
		     mobss(p2)+pstats*sobss(p2),'LineWidth',1);
        hold off
        % Label around the truths; note that Figure 1's tick labels are about
        % the mean (mobss) not the truth (th0), differently here
        set(ah(ind),'xtick',round(rondo*[th0(p1)+vstats*sobss(p1)])/rondo,...
	   'ytick',round(rondo*[th0(p2)+vstats*sobss(p2)])/rondo)

        axis square;
        %grid on
        % To only show lines at +/-2 std and the truth
        hold on
        xlim(th0(p1)+tstats*sobss(p1))
        ylim(th0(p2)+tstats*sobss(p2))
        e2{ind}=plot(repmat(round(rondo*[th0(p1)+...
                    vstats([1 3 5]).*sobss(p1)])./rondo,2,1),...
                    ylim,'-','Color',grey);
        f2{ind}=plot(xlim,...
                    repmat(round(rondo*[th0(p2)+...
                    vstats([1 3 5]).*sobss(p2)])./rondo,2,1),...
                     '-','Color',grey);
        % Color mix
        cmix=[0 0 0]; cmix([p1 p2])=1/2;
        set([p(ind) m(ind)],'MarkerFaceColor',cmix,'MarkerEdgeColor',cmix,'MarkerSize',2)
        % Dull colors
        set([p(ind) m(ind)],'MarkerFaceColor',grey,'MarkerEdgeColor',grey,'MarkerSize',2)

        % Cosmetix
        % Emphasize the OBSERVED means and standard deviations
        set([o1(ind) o2(ind)],'LineWidth',1,'Color',grey(3.5))
        set( m(ind)           ,'MarkerFaceColor',grey(3.5),'MarkerEdgeColor',grey(3.5),'MarkerSize',4)
        top(m(ind),ah(ind))
        bottom(o1(ind),ah(ind))
        bottom(o2(ind),ah(ind))

        % Delete the little cross
        delete([c1(ind) c2(ind)])
        % Delete the other cross
        delete([t1(ind) t2(ind)])
        % Labels
        xl2(ind)=xlabel(flabs{p1});
        ylabel(flabs{p2})

        % Plot pairwise error ellipses
        % https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
        hold on
        % Compute the eigenvectors and eigenvalues of the covariance
        % Think of the Schur complement? How does that relate?
        [V,D]=eig(cov(thhats(:,[p1 p2])));
        % The below is demonstrably not it for blurs=Inf but acceptable for blurs=-1
        % [V,D]=eig(covavhs([p1 p2],[p1 p2])./[sclth0([p1 p2])'*sclth0([p1 p2])]);
        a=sqrt(s)*V*sqrt(D)*[cos(t); sin(t)];
        ep(ind)=plot(a(1,:)+mobss(p1),a(2,:)+mobss(p2));

        % Count the number of estimates outside the confidence interval
        cli=round(sum(inpolygon(thhats(:,p1),thhats(:,p2),...
                                a(1,:)+mobss(p1),a(2,:)+mobss(p2)))/size(thhats,1)*100);
        disp(sprintf('CL asked %g ; received %g ',cl*100,cli))
        
        % And for the calculated covariance, too
        znp=zeros(1,np);
        znp([p1 p2])=1;
        [V,D]=eig(matslice(cvth./(sclth0'*sclth0),znp));
        a11=sqrt(s)*V*sqrt(D)*[cos(t);sin(t)];
        a=sqrt(s)*V*sqrt(D)*[cos(t); sin(t)];
        ec(ind)=plot(a(1,:)+mobss(p1),a(2,:)+mobss(p2));
        % seemax([ah(1) ah(2)],1)
        % seemax([ah(2) ah(3)],2)

        % Dull colors
        set(ep(ind),'LineWidth',1.5,'Color',grey(3.5))
        set(ec(ind),'LineWidth',1,'Color','k')
        % Send these ellipses to the back so the dots show on top
        bottom(ec(ind),ah(ind))
        bottom(ep(ind),ah(ind))
        % Send the grid lines all the way to the back
        for jnd=1:length(e2{ind})
            bottom(e2{ind}(jnd),ah(ind))
            bottom(f2{ind}(jnd),ah(ind))
        end
        if ind==1
            axlm(1,:)=xlim;
            axlm(2,:)=ylim;
        elseif ind==2
            axlm(3,:)=ylim;
        end
        axis square
        xel=xlim; yel=ylim;
        ttt(ind)=text(xel(1)+range(xel)/nrx/2,yel(1)+range(yel)/nry,sprintf('%i%%',cli));
        
        hold off
        longticks(ah)

        if ind==1
            titi=ostitle(ah,params,name,length(thhats(:,1))); movev(titi,-2)
            set(titi,'FontSize',12);
        end
        try
            tt=supertit(ah(1:np),sprintf('%s\n%s%s',sprintf(answs,answ{:}),pm,...
                                         'one sigma uncertainty based on the ensemble'));
            movev(tt,-7)
        end
    end
    % Align the x-axis labels
    for ind=1:3
        xl2(ind).VerticalAlignment='bottom';xl2(ind).Units='centimeters';
    end
    movev(xl2(3),-1);
    xl2(1).Position(2)=xl2(3).Position(2);xl2(2).Position(2)=xl2(3).Position(2);
    moveh(ah(1),-0.02); moveh(ah(3),0.02);

    % Late-breaking
    if np>3; movev(tt,-4); else; movev(tt,-3.5); end
    set(tt,'FontSize',10)
    set(ah,'FontSize',12)
    
    if any(ifinv~=[1 1 1])
        delete(ah(~~ifinv([3 2 1])))
    end
end 

% Output
varns={ah,ha,yl,xl,tl,p,m,o1,o2,ep,ec,axlm,cvth,ttt};
varargout=varns(1:nargout);
