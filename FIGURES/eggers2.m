function eggers2
% EGGERS2
%
% Makes a FIGURE the energy spectral density, the spatial autocovariance and
% simulated fields from isotropic Matern processes. The main routines being
% tested are MATERNOS, MATERNOSY, and SIMULOSL.
%
% Tested on 8.3.0.532 (R2014a) and 8.5.0.197613 (R2015a)
% and 9.0.0.341360 (R2016a) - mileage does vary
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2018
% Last modified by olwalbert-at-princeton.edu, 12/09/2025

% Option for loglog sdf
llsdf=1;
 
% Metric conversion
mfromkm=1000;
% Variance of the "topography" [m^2]
S2=[1^2    2^2    3^2]*mfromkm^2;
% Mean-squared differentiability
NU=[0.5    1     2.0];
% Correlation parameter [m]
RH=[ 250   250    175]*mfromkm;
% The fluctuation scale [m^4]
S0=S2.*pi.*RH.^2/4;

% Adjust manually to where the grid lines overlap the labels in row 2
wout=3000;

% Wavenumber axes, 2D axes flattened to 1D radial
Nx=128; 
dx=30*mfromkm;
% Wavenumber limits [rad/m !]
[k,~,~,dci]=knum2([1 Nx],[1 (Nx-1)*dx]);
k=k(dci(2):end);
% We quote the S(k) on a MUCH finer grid
k=unique([linspace(k(1),k(2),50)...
	  linspace(k(2),k(3),50)...
	  linspace(k(3),k(4),50)...
	  k(4:end)]);
% Spatial axes, for the distance [m]
x=linspace(0,(Nx-1)*dx,100);
% New range for Ky
if llsdf
  newr1=[1e-5 1.75];
  newr2=[-0.1 1.1];
else
  newr1=[-0.1 1.1];
  newr2=[-0.1 1.1];
end
% New ticks for Ky
newy=[0 1/3 2/3 1];

% Major tick wavelengths [km]
xtil=[100 1000 10000];

% Wavenumber limits [rad/km !]
klims=[k(1+1+1) max(k)]*mfromkm;
% Ticks for the log x axes [rad/km !]
exes=10.^[-4 -3 -2 -1];

% Location of labels in wavelength space
xa=1.1*2*pi/0.01;
% Location of labels in wavenumber space
xla=8e-5;
% Location of labels in spatial space
xb=0.975*range(x)/mfromkm;

% Figure axes
clf
fig2print(gcf,'landscape')
ncol=3;
[ah,ha]=krijetem(subnum(length(S2),ncol));

for ind=1:length(S2)
  % Assign the values in question
  th0=[S2(ind) NU(ind) RH(ind)];
  % The normalized spectral density%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(ind))
  Sk= maternos(k,th0,[],2);
  % Conversion is happening 
  if llsdf
    pS(ind)=loglog(k*mfromkm,Sk/S0(ind));
  else
    pS(ind)=semilogx(k*mfromkm,Sk/S0(ind));
  end
  % Keep track of maximum for plotting only
  hold on
  % Do this before the extra axis, the plus one is crucial
  if llsdf
      set(ah(ind),'xlim',klims)
      pgy{ind}=plot(newr1,repmat([exes(3:end) 1],2,1),'Color',grey(9));
  else
    set(ah(ind),'xlim',klims,'ygrid','on')
  end
  % Omit the last one since it's so close to the axis
  pg{ind}=plot([exes(1:end-1) ; exes(1:end-1)],newr1,'Color',grey(9));
  % These are the percentages
  pers=[25 50 75];
  % Remove the ones where the labels will be
  if llsdf
      yremr=[0.03 0.3];
      try
          yreml=[0.5 1.5].*10.^(-3+0.5.*[1:length(pers)]');
      catch
          yreml=repmat([0.5 1.5],3,1).*repmat(10.^(-2.9+0.5*[1:length(pers)]'),1,2);
      end
  else
      yremr=[0.41 0.99];
      try 
          yreml=[0.01 0.19]+([1;2;3]-1)*0.2;
      catch
          yreml=repmat([0.01 0.19],3,1)+repmat([1:length(pers)]'-1,1,2)*0.2;
      end
  end
  % White vertical for the parameter cover
  pgg{ind}=plot([exes(end-1) ; exes(end-1)],yremr,'w-');
  if llsdf
    % White horizontal
    pgyg(1)=plot([exes(1)*0.75 exes(2)*1.5],[exes(end-1) ; exes(end-1)],'w-');
    pgyg(2)=plot([exes(3)*0.75 exes(end)],[exes(end) ; exes(end)],'w-');
  end
  % Plot some Matern percentages!
  for ond=1:length(pers)
    Kprc(ond)=maternprc(th0,pers(ond));
    hold on
    ppr{ind}(ond)=plot([Kprc(ond) Kprc(ond)]*mfromkm,newr1,'k-');
    % White vertical
    pww{ind}(ond)=plot([exes(1) exes(1)],yreml(ond,:),'w-',...
                       'linew',2);
    % White vertical
    pwx{ind}(ond)=plot([exes(2) exes(2)],yreml(ond,:),'w-',...
                       'linew',1);
    % Annotate the percentiles
    if llsdf
      yla=10^(-3+0.5*ond);
    else
      yla=0.1-0.035+(ond-1)*0.2;
    end
    al(ond,ind)=text(xla,yla,...
                     sprintf('%s_{%i} = %i km','\lambda',...
                             pers(ond),round(2*pi/Kprc(ond)/mfromkm)));
  end
  if llsdf
    ylim([1e-3 newr1(2)])
  else
    ylim(newr1)
  end
  hold off

  if ind==1
    yl(1)=ylabel('normalized sdf S(k)/S(0)','FontSize',9);
  end
  xl(ind)=xlabel('wavenumber (rad/km)','FontSize',9);
  
  % Put on the wavelength axis
  set(ah(ind),'xtick',exes)
  [axx(ind),xl(ind),yl(ind)]=...
      xtraxis(ah(ind),xtil,xtil,...
	      sprintf('wavelength (km)'));
  set(axx(ind),'xdir','rev','xlim',2*pi./klims([2 1]))
  % Minor visual adjustment that I was not able to correct anywhere else
  shrink(axx(ind),1,0.99); movev(axx(ind),-0.001)
  % Legends on the extra axis
  axes(axx(ind))
  if llsdf
    ya=[0.2 0.1 0.05];
  else
    ya=[0.9 0.7 0.5];
  end
  aa(1,ind)=text(xa,ya(1),sprintf('%s = %3.1f km','\sigma',sqrt(S2(ind))/mfromkm));
  aa(2,ind)=text(xa,ya(2),sprintf('%s = %3.2f %s','\nu',NU(ind),'  '));
  aa(3,ind)=text(xa,ya(3),sprintf('%s = %i km','\rho',RH(ind) /mfromkm));

  % Reorder
  for ond=1:length(ppr)
      top(ppr(ond),ah(ind))
  end
  top(pS,ah(ind))
  
  % The correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(ind+ncol))
  Kx=maternosy(x,th0);
  pK(ind)=plot(x/mfromkm,Kx/S2(ind));
  set(ah(ind+ncol),'xlim',[x(1) max(x)]/mfromkm,...
                   'ytick',newy,'yticklabel',round(newy*10)/10,...
                   'ygrid','on','xgrid','on')
  hold on
  % Remove the ones where the labels will be
  pgw{ind}=plot([wout wout],[0.81 0.99],'w-','linew',2);
  % At this point the correlations should have died down by one third
  pR(ind)=plot(pi*[1 1]*RH(ind)/mfromkm,newr2,'b');
  % At this point the limit of the simulation has been reached
  % pL(ind)=plot([1 1]*(Nx-1)*dx/mfromkm,newr,'k');
  hold off
  bb(1,ind)=text(xb,0.9,sprintf('%s = %i km',...
                                '\pi\rho',round(RH(ind)*pi/mfromkm)));
  
  if ind==1
    yl(2)=ylabel('correlation C(r)/\sigma^2','FontSize',9);
    yl(2)=ylabel('normalized cov C(r)/\sigma^2','FontSize',9);
  end
  xl(ind+ncol)=xlabel('distance (km)');
  
  % The fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(ind+2*ncol))
  % Sampling steps
  params.dydx=[dx dx];
  % Sizes of the fields
  params.NyNx=[Nx Nx];
  % Spectral blurring using BLUROS or BLUROSY
  params.blurs=Inf;
  % Space decorrelation using QUARTER
  % params.quart=1; disp('Quartering')
  % Isotropic filtering
  % params.kiso=pi/sqrt(prod(params.dydx)); disp('Isotropic filtering')
  % params.kiso=NaN; disp('(Rather) no filtering')
  % Perform the actual simulation
  Hx=simulosl([S2(ind) NU(ind) RH(ind)],params);
  % Plot it up to within twice the actual or sample variance
  Hxm=mean(Hx);
  Hx=Hx-Hxm;
  % Knowing this isn't normal!! Unless you nix the correlation
  display(sprintf('sample std %3.1f km population std %3.1f km',...
          std(Hx)/mfromkm,sqrt(S2(ind))/mfromkm));
  % What  you think it should be... probably rather than what it is?
  sox=sqrt(S2(ind)); disp('Using the ACTUAL variance')
  sox=std(Hx); disp('(Rather) using the SAMPLE variance')
  sax=[-2 2]*sox;
  ahh(ind)=imagefnan([0 0],[(Nx-1)*dx (Nx-1)*dx]/mfromkm,v2s(Hx),[],sax);
  axis image
  % Plot a circle somewhere that shows the decorrelation region
  hold on
  hc(ind,1:2)=circ(pi*RH(ind)/mfromkm,[],[1000 1000]); delete(hc(ind,2))
  hold off
  tl(ind)=title(sprintf('  m = %3.1f km ; s = %3.1f km',...
                        Hxm/mfromkm,std(Hx)/mfromkm));
  xl(ind+2*ncol)=xlabel('easting (km)');
  if ind==1
    yl(3)=ylabel('northing (km)','FontSize',9);
  end
  % Reasonable color axis in terms of standard deviations
  gp=getpos(ah(ind+2*ncol));
  [cb(ind),xcb(ind)]=addcb([gp(1)+gp(3)*0.95 gp(2)*0.99 gp(3)/12 gp(4)],...
                           sax,sax,'kelicol',sox);
  hold on
  pgk(ind)=plot(get(cb(ind),'xlim'),[0.5 0.5],'k');
  hold off
  axes(cb(ind))
  xxb=2*indeks(get(cb(ind),'xlim'),2);
  yyb=get(cb(ind),'ytick');
  sl(1)=text(xxb,yyb(1),sprintf('%s','-2s'));
  sl(2)=text(xxb,yyb(2),sprintf('%s','-1s'));
  sl(3)=text(xxb,yyb(4),sprintf('%s','+1s'));
  sl(4)=text(xxb,yyb(5),sprintf('%s','+2s'));
  sl(5)=text(xxb,yyb(3),sprintf('%s','m'));
  nolabels(cb(ind),2)
  delete(xcb)  
end

% Cosmetics
longticks([ah(1:2*ncol) axx],1)
longticks(ah([2*ncol+1:3*ncol]),0.75)
set(aa,'HorizontalAlignment','l')
set(al,'HorizontalAlignment','l','FontSize',9)
set([tl xl yl],'FontSize',9)
set(bb,'HorizontalAlignment','r','Color','b')
set(hc(:,1),'color','b')
set(ah([1:length(S2)]+3),'ylim',newr2)
nolabels(ha(4:end),2)
set(pS,'linew',2,'color','k')
set(pK,'linew',2,'color','b')
moveh([ah(2*ncol+1:3*ncol) cb(1:end)],-0.0175)
moveh([ah(2*ncol+2:3*ncol) cb(2:end)],-0.0050)
moveh(tl,250)

% Printo
figna=figdisp([],[],[],1);
keyboard

system(sprintf('epstopdf %s.eps',figna));
