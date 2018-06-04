function eggers3(domind,fldind)
% EGGERS3(domind,fldind)
%
% A whole suite of runs for the Venus topography!
%
% EXAMPLE
%
% parfor index=1:77; try; eggers3(index); end; end
% for index=1:77; try; eggers3(index); end; end
% [parfor gets in trouble with the graphics]
%
% Last modified by fjsimons-at-alum.mit.edu, 07/14/2015

% needs some work 

% Plot patch?
img=0;

% Plot residuals? Turn on for now
chi=0;

% Simulate some relateds?
sos=0;

% Plot covariance in space and in the spectrum?
kov=1;

% Estimate uncertainties? Turn off for now
unc=0;

% Sets the domain index
defval('domind',58);
% Sets the field in question [03 for 360 topography]
defval('fldind',3);

% Make the variable name and the save file to look for
varibal=sprintf('V%4.4i_%2.2i',domind,fldind);
fname=sprintf('%s.mat',varibal);

% Load the data, take a look
topo=modload(fullfile(getenv('IFILES'),'VENUS','plmData',...
                      'plmVenus_D-5.mat'),varibal);
% Reference one-degree equatorial lontitude distance
DegDis=2*pi*fralmanac(topo.id.fralRef,'Venus')/360;

% Harmonize the data names
Hx=topo.dataP.dp(:);

% Prepare for color bar and labeling
cm=2;
cax=[-cm cm]*std(Hx);
cxl=[-cm:cm]*std(Hx);

% Check the spacings etc
diferm(abs(diff(topo.geo.lonrDx))/[topo.params.NyNx(2)-1],topo.geo.DxDy)
diferm(abs(diff(topo.geo.latrDx))/[topo.params.NyNx(1)-1],topo.geo.DxDy)
c11=[topo.geo.lonrDx(1) topo.geo.latrDx(2)];
cmn=[topo.geo.lonrDx(2) topo.geo.latrDx(1)];

% PERFORM THE ESTIMATION!
if exist(fullfile('.',fname))~=2 
  % Complete the parameters, make blurring 3 for safety
  topo.params.blurs=3;
  % No need for the quartering here, by the way
  topo.params.quart=0;
  % Put in a good choice for kiso as the Nyquist wavenumber
  topo.params.kiso=pi./max(topo.params.dydx);
  %topo.params.kiso=NaN;

  if ~isnan(topo.params.kiso)
    disp('May need to widen the default bounds on nu inside MLEOSL')
    disp('May need to switch to unconstrained search inside MLEOSL')
  end

  % Estimates the parameters via maximum-likelihood... using a reasonable
  % but randomized initial guess! Stick to constrained inversions for now
  % but widen the bounds, especially for nu
  algo='con';
  [thhat,covh,logli,thini,scl,p,e,o,gr,hs,Hk,k,ops,bds]=...
      mleosl(Hx,[],topo.params,algo);

  % Here we compute and write out the moments of the Xk
  [L,~,momx,vr]=logliosl(thhat,p,Hk,k,scl);
  [pf,pv,pr,af]=normtest(momx(3),1,vr,0.05);
  disp(sprintf('NORMTEST %i %5.3f %i',pf,pv,round(pr)))
  if pf==0; stp='accept'; else stp='reject'; end

  % Scale the parameters back to physical units
  thhat=thhat.*scl;

  % What do we expect for var(Hx), biased as it is with this rho?
  b=varbias(thhat,topo.params);

  % Take a sanity-check look at the variance for comparison
  disp(sprintf('var(Hx) = %i ; s^2 = %i ; debiased %i',...
               round(var(Hx)),round(thhat(1)),round(thhat(1)-b)))

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SAVE ALL OF THE IMPORTANT STUFF (EXCEPT THE DATA AGAIN)
  save(varibal,'thhat','covh','logli','thini','p','e','o','gr','hs','b',...
       'momx','vr','stp','pf','pv','af','Hk','varibal','DegDis','ops','bds')
else
  load(fname)
end

% Pull out the parameters for easy reference
s2=thhat(1);
nu=thhat(2);
rh=thhat(3);

% Plot data patch
if img==1
  % Makes a plot of the data in question
  clf; fig2print(gcf,'portrait')
  ah=gca;
  [hp,cp]=imagefnan(c11,cmn,topo.dataP.dp-mean(topo.dataP.dp(:)),...
               'kelicol',cax); hold on
  xl(1)=xlabel(sprintf('apparent longitude [%s]',str2mat(176)));
  yl(1)=ylabel(sprintf('apparent latitude [%s]',str2mat(176)));
  longticks(ah,2)
  % Plots the region on top of it
  hr=twoplot(topo.geo.XYr360,'LineW',2,'Color','k');
  hold off

  % Color bar
  [cb,xb]=addcb('vert',cax,cax,'kelicol',diff(cxl));
  set(cb,'ytickl',round(cxl),'yaxisl','r');
  set(xb,'string',sprintf('demeaned topography of %s [m] (to within %i std)',...
                          nounder(varibal),cm))
  longticks(cb,2)
  drawnow

  % Make a title 
  axes(ah)
  t=title(sprintf(...
      '%s  |  %s = %i   %s = %4.2f   %s = %i km  |  %s = %i km  |  %s  pval = %4.2f',...
      nounder(varibal),...
      '\sigma',round(sqrt(s2)),'\nu',nu,'\rho',round(rh/1000),...
      '\lambda',round(2*pi/p.kiso/1000),stp,pv));
  set(t,'FontS',12)
  movev(t,range(topo.geo.latrDx)/50)

  % Save in high-quality format
  figna=figdisp([],varibal,[],2);
end

% Plot residual behavior
if chi==1
  clf
  % We may have had no wavenumber restriction; report WITH restriction
   p.kiso=pi./max(p.dydx);
   mlechipsdosl(Hk,thhat,p,...
                sprintf('%s  |  %s = %i m  %s = %4.2f  %s = %i km',nounder(varibal),...
                        '\sigma',round(sqrt(s2)),'\nu',nu,'\rho',round(rh/1000)))
   figna=figdisp([],sprintf('%s_chi',varibal),[],2);
end

% How about some simulations with similar parameters
if sos==1
  clf
  [ah,ha,H]=krijetem(subnum(2,2));
  axes(ah(1))
  h=imagefnan(c11,cmn,topo.dataP.dp-mean(topo.dataP.dp(:)),...
                       'kelicol',cax); hold on
  t(1)=title(sprintf('%s',nounder(varibal)));

  [cb,xb]=addcb('vert',cax,cax,'kelicol',diff(cxl));
  set(cb,'ytickl',round(cxl),'yaxisl','r');
  set(xb,'string',nounder(varibal))

  % We may have had blurring; simulate without blurring? Or do we 
  % I go back and forth on this one but always make sure that kiso is in there!
  p.blurs=0;
  for ind=2:4
    Hx=simulosl(thhat,p);
    axes(ah(ind))
    h(ind)=imagefnan(c11,cmn,reshape(Hx,p.NyNx)-mean(Hx(:)),...
                       'kelicol',cax);
    [cb(ind),xb(ind)]=addcb('vert',cax,cax,'kelicol',diff(cxl));
    set(cb(ind),'ytickl',round(cxl),'yaxisl','r');
    set(xb(ind),'string',sprintf('simulation %i',ind-1))
  end
  hold off
  longticks(ah)
  axes(ah(2))
  t(2)=title(sprintf(...
      '%s = %i m  %s = %4.2f   %s = %i km',...
      '\sigma',round(sqrt(s2)),'\nu',nu,'\rho',round(rh/1000)));

  % Cosmetics
  nolabels(ah(1:2),1)
  nolabels(ha(3:4),2)
  serre(H,2/3,'across')
  serre(H',2/3,'down')
  delete(cb([1 2 3]))
  delete(xb(4))

  figna=figdisp([],sprintf('%s_sim',varibal),[],2);
end

if kov==1
  % Inspired by EGGERS2
  % Wavenumber axes, 2D axes flattened to 1D radial
  Nx=sqrt(prod(p.NyNx)); 
  dx=sqrt(prod(p.dydx)); ;
  % Wavenumber limits [rad/m !]
  [k,~,~,dci]=knum2([1 Nx],[1 (Nx-1)*dx]);
  k=k(dci(2):end); kor=k;
  % We quote the S(k) on a MUCH finer grid
  k=unique([linspace(1/Nx/dx,k(1),50)...
            linspace(k(1),k(2),50)...
	    linspace(k(2),k(3),50)...
	    linspace(k(3),k(4),50)...
	    k(4:end)]);
  % Spatial axes, for the distance [m]
  exor=linspace(0,(Nx-1)*dx,Nx);
  xpnd=3;
  x=linspace(0,xpnd*(Nx-1)*dx,100);
  % Adjust manually to where the grid lines overlap the labels in row 2
  wout=(xpnd-1)*(Nx-1)*dx/1000;
  % New range for Ky
  newr=[-0.1 1.1];
  % New ticks for Ky
  newy=[0 1/3 2/3 1];

  % Major tick wavelengths [km]
  xtil=[1e3 1e4 1e5];
  xtill={'1e3' '1e4' '1e5'};
  
  % Wavenumber limits [rad/km !]
  klims=[k(1+1+1) max(k)]*1e3;
  % Ticks for the log x axes [rad/km !]
  exes=10.^[-4 -3 -2];
  % Spatial limits [km !]
  xlims=[x(1) x(end)]/1e3;
  
  % Ticks for the spatial x axis [km !]
  Kexes=round([0 cumsum(repmat((Nx-1)*dx,1,xpnd))]/1000);

  % Location of labels in wavelength space
  xa=1.1*2*pi/0.01*5;
  % Location of labels in spatial space
  xb=0.975*range(x)/1e3;
  
  % Calculate the covariance
  Sk= maternos(k,thhat,2);
  Skor=maternos(kor,thhat,2);

  Kxor=maternosy(exor,thhat,2);
  Kx=maternosy(x,thhat,2);
  
  clf; fig2print(gcf,'portrait')
  [ah,ha]=krijetem(subnum(1,3));
  % Plot the spectral covariance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  % Conversion is happening 
  pSk=semilogx(k*1000,Sk/[thhat(1).*pi.*thhat(3).^2/4]); hold on
  pSor=semilogx(kor*1000,Skor/[thhat(1).*pi.*thhat(3).^2/4],'o');
  % Remove the last one since it's so close to the axis
  pg{1}=plot([exes(1:end-1) ; exes(1:end-1)],newr,':');
  % Remove the ones where the labels will be
  pgg{1}=plot([exes(end-1) ; exes(end-1)],[0.41 0.99],'w-');
  hold off
  yl(1)=ylabel('normalized sdf S(k)/S(0)');
  xl(1)=xlabel('wavenumber (rad/km)');

  set(ah(1),'xlim',klims,'ygrid','on','xtick',exes)
  [axx(1),xl(4),yl(4)]=...
      xtraxis(ah(1),xtil,xtill,...
	      sprintf('wavelength (km)'));
  set(axx(1),'xdir','rev','xlim',2*pi./klims([2 1]))

  % Legends on the extra axis
  axes(axx(1))
  aa(1)=text(xa,0.975,sprintf('%s = %i m','\sigma',round(sqrt(thhat(1)))));
  aa(2)=text(xa,0.775,sprintf('%s = %4.2f %s','\nu',thhat(2),'  '));
  aa(3)=text(xa,0.575,sprintf('%s = %i km','\rho',round(thhat(3)/1e3)));
  
  % Plot the spatial covariance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  pK=plot(x/1000,Kx/thhat(1)); hold on
  pKor=plot(exor/1000,Kxor/thhat(1)); 
  set(ah(2),'xlim',xlims,'ygrid','on')
  % Remove the first and last one since they are so close to the axis
  pw{1}=plot([Kexes(2:end-1) ; Kexes(2:end-1)],newr,':');
  % Remove the ones where the labels will be
  pgw{1}=plot([wout wout],[0.81 0.99],'w-','linew',2);
  % At this point the correlations should have died down by one third
  pR=plot(pi*[1 1]*thhat(3)/1000,newr);
  bb(1)=text(xb,0.9,sprintf('%s = %i km',...
                                '\pi\rho',round(thhat(3)*pi/1000)));

  yl(2)=ylabel('correlation C(y)/\sigma^2');
  xl(2)=xlabel('distance (km)');

  % Plot the field again, or one realization thereof %%%%%%%%%%%%%%%%%%%%
  axes(ah(3))
  h=imagefnan(c11,cmn,topo.dataP.dp-mean(topo.dataP.dp(:)),...
                       'kelicol',cax); hold on
  % Center this panel, knowing that its dimensions do change
  layout(ah(3),0.5,'m','y')
  xl(3)=xlabel(sprintf('apparent longitude [%s]',str2mat(176)));
  yl(3)=ylabel(sprintf('apparent latitude [%s]',str2mat(176)));
  t=title(sprintf('%s',nounder(varibal)));
  
  % Cosmetics 
  set(aa,'HorizontalA','l')
  set(bb,'HorizontalA','r','Color','b')
  set([ah(1) axx(1)],'position',...
            [getpos(ah(1),1) 0.375 getpos(ah(1),3) 0.25])
  set(ah(2),'position',...
            [getpos(ah(2),1) 0.375 getpos(ah(2),3) 0.25])
  set(ah(1:2),'ylim',newr,'ytick',newy,'ytickl',round(newy*10)/10)
  set(ah(2),'xtick',Kexes)
  % Show the wavenumbers that we actually have 
  set(pSor,'color','k','marker','o','markers',3,'markere','k','markerf','k')
  % Show the maximum range that we could ever observe
  set(pKor,'color','b','LineW',2)
  set(pSk,'linew',1,'color',grey)
  set(pK,'linew',1,'color',grey)

  longticks([ah axx])
  set(findobj('fontsize',10),'fontsize',8)
  set([xl yl],'fontsize',8)

  figna=figdisp([],sprintf('%s_kov',varibal),[],2);
  % system(sprintf('xpdf %s.pdf',figna))
end

if unc==1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ESTIMATE UNCERTAINTIES
  
  % Get an idea of the uncertainty by multiply running the whole thing again
  % for made-up data on the basis of the recovered parameters
  spmd
    try
      % If it "hits the wall" it won't hit the spot
      mleosl('demo1',8,thhat,topo.params)
    catch
      % In which case we need to keep going (or close out a run?)
    end
  end
  % Rename all those files
  system(sprintf('osrename %s %s_%s',date,date,varibal));
  % And then work with these guys to come up with empirical covariances
  thhats=load(sprintf('mleosl_thhat_%s_%s',date,varibal))
  % Derive a new covariance
  covhs=cov(thhats);
  % Compare with the other covariances
end
