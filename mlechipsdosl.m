function varargout=mlechipsdosl(Hk,thhat,scl,params,stit,ah)
% [a,mag,ah,ah2,cb,ch,spt]=MLECHIPSDOSL(Hk,thhat,scl,params,stit,ah)
%
% Makes a four-panel plot of the quadratic residuals and their
% interpretation for a Matern likelihood model of a single-field
% Fourier-domain data patch, as well as plots of the actual power spectral
% density of that data patch and its predicted power spectral density based
% on the model. These residuals should have a chi-squared distribution, so
% the produced figure is useful in determining how well the Matern
% likelihood model fits the data.
% 
% The figure produced consists of four panels:  
%    (TL) Histogram and qq-plot of the chi-squared residuals
%    (TR) Predicted 2D power spectral density based on "thhat"
%    (BL) 2D chi-squared residuals
%    (BR) "True" power spectral density with contours from (TR) overlain
%
% INPUT:
%
% Hk         The Fourier-domain data, e.g. from SIMULOSL
% thhat      The evaluated scaled Matern parameter vector
% scl        The scale factors
% params     The structure with the fixed parameters from the experiment
% stit       Title for the overall figure [defaulted]
% ah         A quartet of axis handles [defaulted]
%
% OUTPUT:
%
% a           1 or 0 whether the model was accepted or not
% mag         The "magic" parameter which I've come to call s_X^2
% ah,ah1      Main and secondary axis (panels 1-4) handles
% cb          Colorbar (panel 3) handle
% ch          Contour (panels 2 and 4) handles
% spt         Supertitle handle
%
% SEE ALSO:
%
% EGGERS8, MATERNOS, LOGLIOSL, SIMULOSL, QQPLOT
%
% NOTE: 
%
% Maybe should integrate MLECHIPLOS into this one. 
%
% Last modified by gleggers-at-princeton.edu, 04/17/2014
% Last modified by fjsimons-at-alum.mit.edu, 08/24/2017

% Some defaults
defval('stit','Chi-squared residuals')
defval('ah',krijetem(subnum(2,2)))

% The following is as in MLECHIPLOS

% So now you have found a solution and you evaluate it at thhat
[~,~,~,k,~,Sb,Lb]=simulosl(thhat.*scl,params,0);
% Get just a little more information
[kk,~,~,kx,ky]=knums(params); 

% Degrees of freedom for 2X system (assuming a single field)
df=2;

% Get the quadratic residuals, any of a number of ways

% Multiply to obtain a variable which should follow the rules 
%Zk=[Hk(:)./Lb(:)];
%Xk0=hformos(1,Zk,[1 0 1]);
% Same thing
Xk=abs(Hk(:)).^2./Sb(:);

% Calculate residuals, removing values at the zero wavenumber and
% wavenumbers above the minimum wavenumber "params.kiso"
Xk(k==0)=NaN;
Xk(k>params.kiso)=NaN;

% And this should be the same thing again, except how it treats k=0
[Lbar,~,~,momx,~,Xk1]=logliosl(k,thhat,scl,params,Hk,1);
Xk1=-Xk1-log(Sb(~~k));
% The oldest way, using a since retired function
% Xkk1=-Lkosl(k,thhat.*scl,params,Hk)-log(Sb);
% difer(Xkk1(~~k)-Xk1,9,[],NaN)

% Check we're doing the same thing to tolerance, depending on whether
% some prior codes put a NaN at zero wavenumber or got rid of the
% zero-wavenumber values altogether; these checks will be removed
%difer(Xk(~isnan(Xk0))-Xk0(~isnan(Xk0)),9,[],NaN)
%difer(Xk(~~k)-Xk1,9,[],NaN)

% Labeling  thing
varibal='X';
xstr2v=sprintf('quadratic residual 2%s',varibal);

% Evaluate the likelihood
% disp(sprintf('The loglihood is %8.3f',Lbar))

%% TAKE CARE OF A FEW THINGS UPFRONT
% Bounds of X residuals to show (functions as color axis in panel 3)
boundsX0=[0 3*df];
% Bounds of 2X residuals to how (functions as axis bounds in panel 1)
bounds2X=2*boundsX0;
% Color axis for the power spectral densities (panels 2 and 4)
caxPSD=[-5 0];
% Contours to be plotted on the power spectral density plots
conPSD=[-5:-1];
% Get order of magnitude of last wavenumber for scaling
om=round(log10(max(k(:))));

% Set up the wavenumber/wavelength axes labels for panels 2-4
xstr1=sprintf('x wavenumber (rad/m) %s10^{%i}','\times',om);
ystr1=sprintf('y wavenumber (rad/m) %s10^{%i}','\times',om);
xstr2='x wavelength (km)';
ystr2='y wavelength (km)';

% Prepare the overall figure
fh=gcf;
fig2print(fh,'landscape')
set(fh,'PaperPositionMode','auto')

% Select non-ridiculous values
if any(isinf(Xk)); error('Adapt as in MLECHIPLOS and reconcile'); end
allg=~isnan(Xk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PANEL 1: TWICE CHI-SQUARED RESIDUALS AND QQ-PLOT [Top-left]
axes(ah(1))

% Get binning info for the twice chi-squared residuals histogram
binos=5*round(log(length(Xk(allg))));
% Rather this, to hit the grid lines
binWidth=df/4;
binos=(binWidth/2)*[1:2:bounds2X(2)*2/binWidth+1];
[bdens,c]=hist(2*Xk(allg),binos);
% Plot the histogram as a bar graph
bdens=bdens/indeks(diff(c),1)/length(Xk(allg));
bb=bar(c,bdens,1);
% Each of the below should be df/2
%disp(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
%		   varibal,nanmean(Xk(allg)),...
%		   varibal,nanvar(Xk(allg))));

set(bb,'FaceC',grey)
hold on
% Plot the ideal chi-squared distribution
refs=linspace(0,max(bounds2X),100);
plot(refs,chi2pdf(refs,df),'Linew',1.5,'Color','k')
hold off

% Labeling and cosmetic adjustments
xlim(bounds2X)
xl1(1)=xlabel(xstr2v);
ylim([0 max(bdens)*1.05])
yl1(1)=ylabel('probability density');

% Prepare an overlay axis for the quantile-quantile plot
ah2(1)=laxis(ah(1),0,0);
axes(ah2(1))

% Obtain qq-plot data
h=qqplot(2*Xk,ProbDistUnivParam('gamma',[df/2 2]));
hx=get(h,'Xdata'); hx=hx{1};
hy=get(h,'ydata'); hy=hy{1};
delete(h)

% Make the qq-plot
qq0=plot(bounds2X,bounds2X,'k'); hold on
qq=plot(hx,hy,'LineS','none','Marker','o','MarkerF','r',...
	'MarkerE','r','MarkerS',2);

% More cosmetic adjustments
set(ah(1),'box','off')
set(ah2(1),'xlim',get(ah(1),'xlim'),'yaxisl','r','box','off',...
	   'xaxisl','t','color','none','ylim',get(ah(1),'xlim'))

% Add labels and tickmarks
ylr(1)=ylabel('quantile-quantile prediction');
tkmks=bounds2X(1):2*df:bounds2X(2);
set([ah(1) ah2(1)],'xtick',tkmks)
set(ah2(1),'ytick',tkmks)

% Add gridlines  % (c) Jos van der Geest
try
  gh=gridxy([4 8],[4 8],'LineStyle',':');
end

% Information is power
tbstr{1}=sprintf('%5.2f%%',100*sum(binWidth*bdens(1:end-1)));
tbstr{2}=sprintf('max(2X) = %5.2f',max(2*Xk));
tbstr{3}=sprintf('m(X) = %5.2f',nanmean(Xk));
tbstr{4}=sprintf('v(X) = %5.2f',nanvar(Xk));

% Calculate the mean of (X-df/2)^2 so it can be passed as output
magx=nanmean((Xk(allg)-df/2).^2);
%tbstr{5}=sprintf('mean([X-%d]^2) = %5.2f',df/2,magx);
tbstr{5}=sprintf('s_X^2 = %5.2f',magx);

% Do the test whether you accept this as a good fit, or not
vr=8/sum(allg);
[a,b,c,d]=normtest(magx,1,vr,0.05);
% disp(sprintf('NORMTEST %i %5.3f %i',a,b,round(c)))
if a==0; stp='accept'; else stp='reject'; end
tbstr{6}=sprintf(' %s at %4.2f',stp,d);
tbstr{7}=sprintf('p = %5.2f',b);

% Give the x- and y-positions of the textbox
tbx=repmat(bounds2X(2)-bounds2X(2)/30,1,length(tbstr));
tby=[7.5 6.35-[0 1 2 3.5 4.5 5.5]];

% Make the textbox(es) with unbordered fillboxes around them
for i=1:length(tbstr)
  tb(i)=text(tbx(i),tby(i),tbstr{i},'HorizontalA','Right');
  gg=ext2lrtb(tb(i),1.1,1.1); delete(tb(i)); hold on
  fb(i)=fillbox(gg+[-0.1 0.08 0 0]);
  tb(i)=text(tbx(i),tby(i),tbstr{i},'HorizontalA','Right');
  hold off
  set(fb(i),'EdgeC','w','FaceC','w')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PANEL 2: PREDICTED 2D POWER SPECTRUM [Top-right]
axes(ah(2))

% The top-left and bottom-right wavenumbers, scaled
c11=[kx(1) ky(end)]/10^om;
cmn=[kx(end) ky(1)]/10^om;

% Remove the zero wavenumber value (to avoid issues with the
% colorscale) and those at wavenumbers above the isotropic cutoff
Sb(k==0)=NaN;
Sb(k>params.kiso)=NaN;

Sb=reshape(Sb,params.NyNx);

% Scale to the largest predicted value and plot the power spectrum
sclSb=log10(Sb./max(Sb(:)));
imagefnan(c11,cmn,sclSb,'jet',caxPSD,[],[],0);

% Label the wavenumber axes
xl1(2)=xlabel(xstr1);
yl1(2)=ylabel(ystr1);

% For the wavelength axis, get tickmark values on the wavenumber axis
xtk=get(ah(2),'xtick')*10^om;
ytk=get(ah(2),'ytick')*10^om;

% Convert to wavelengths (in m), recognizing the zero wavenumber
xtkl=2*pi./xtk/1000;
xtkl(isinf(xtkl))=[params.NyNx(2)-1]*params.dydx(2)/1000;
ytkl=2*pi./ytk/1000;
ytkl(isinf(ytkl))=[params.NyNx(1)-1]*params.dydx(1)/1000;

% Create and label the wavelength axis
[ah2(2),xl2(2),yl2(2)]=xtraxis(ah(2),xtk/10^om,round(xtkl),xstr2,...
			       ytk/10^om,round(ytkl),ystr2);

% Return to the main axis and prepare to plot contours
axes(ah(2)); hold on

% Calculate and plot contours
[~,ch(1)]=contour(kx/10^om,fliplr(ky/10^om),sclSb,conPSD,'LineW',2);
caxis(caxPSD)

% Option to set the contours to black for visibility
%set(hb,'color','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% PANEL 3: 2D X RESIDUALS [Bottom-left]
axes(ah(3))

% Plot the X residuals
imagefnan(c11,cmn,reshape(Xk,params.NyNx),'gray',boundsX0,[],1,0);

% The BOXTEX messes up the neat flat shading from imagefnan... or
% does it?
if ~isnan(params.kiso)
  % Place a box giving the wavelength of the isotropic wavenumber cutoff
  [~,zy]=boxtex('ll',ah(3),sprintf('%2.0f km',2*pi./params.kiso/1000),...
		12,[],[],1.1);
  set(zy,'fonts',get(zy,'fonts')-1)
end

% Label the wavenumber axes
xl1(3)=xlabel(xstr1);
yl1(3)=ylabel(ystr1);

% Create and label the wavelength axis
ah2(3)=xtraxis(ah(3),xtk/10^om,round(xtkl),xstr2,...
	       ytk/10^om,round(ytkl),ystr2);

% Remove the right y-label to avoid clutter
delete(get(ah2(3),'ylabel'))

axes(ah(3))

% Add a colorbar
[cb,xcb]=addcb('vert',boundsX0,boundsX0,'gray',df,1);
set(xcb,'string','quadratic residual X')

% Reposition the main axis
set(ah(3),'position',[getpos(ah(3),1) getpos(ah(3),2) ...
		    getpos(ah(2),3) getpos(ah(2),4)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PANEL 4: 2D TRUE POWER SPECTRUM [Bottom-right]
axes(ah(4))

% Calculate the periodogram of the edata
Sk=abs(reshape(Hk(:),params.NyNx).^2);

% Remove the zero wavenumber value (to avoid issues with the
% colorscale) and those at wavenubmers above the isotropic cutoff
Sk(k==0)=NaN;
Sk(k>params.kiso)=NaN;

% Scale to the largest predicted value and plot the power spectrum
sclSk=log10(Sk./max(Sb(:)));
imagefnan(c11,cmn,sclSk,'jet',caxPSD,[],[],0);

% Label the wavenumber axes
xl1(4)=xlabel(xstr1);
yl1(4)=ylabel(ystr1);

% Create and label the wavelength axis
ah2(4)=xtraxis(ah(4),xtk/10^om,round(xtkl),xstr2,...
	       ytk/10^om,round(ytkl),ystr2);

% Return to the main axis and prepare to plot contours
axes(ah(4)); hold on

% Place contours from the predicted power spectrum onto this true one
[~,ch(2)]=contour(kx/10^om,fliplr(ky/10^om),sclSb,conPSD,'LineW',2);
caxis(caxPSD)

% Option to set the contours to black for visibility
%set(hb,'color','k')

%% FINAL COSMETIC ADJUSTMENTS

% Adjust tickmarks
longticks([ah ah2 cb])

% Adjust the colorbar in panel 3
axes(cb)
set(cb,'YAxisLocation','right')
moveh(cb,0.06)
shrink(cb,1.3,1.27)

% Adjust the main axes
set(ah(2:4),'Box','On')

% Give the overall figure a title
axes(ah2(1))
%spt=title(ah2(1),stit);
%movev(spt,-.05)
spt=text(df,bounds2X(2)-df+1/2,stit);

movev([ah ah2 cb],-.02)

% Set figure background color to white
set(gcf,'color','W','InvertH','off')

% Collect output
vars={a,magx,ah,ah2,cb,ch,spt};
varargout=vars(1:nargout);
