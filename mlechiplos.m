function varargout=mlechiplos(witsj,Hk,thhat,scl,params,ah,pertur,th0,covX,E,v)
% [cb,xl,t,tt]=MLECHIPLOS(witsj,Hk,thhat,scl,params,ah,pertur)
% [cb,xl,t,tt]=MLECHIPLOS(witsj,Hk,thhat,scl,params,ah,pertur,th0,covX,E,v)
%
% Makes a THREE-panel plot of the quadratic residuals and their
% interpretation for a likelihood model evaluated at its estimate. Plotted
% are a histogram with the theoretical distribution overlain, a
% quantile-quantile plot against the theoretical distribution, and a
% two-dimensional spectral plot to check for patterns. Calculates metrics
% and performs a statistical test for model suitability.
%
% INPUT:
%
% witsj    What program do we actually run?
%          1 SIMULOS - old-school uncorrelated loads
%          2 SIMULROS - new school correlated as in Simons & Olhede (2013)
%          3 SIMULROS0 - new school uncorrelated as in Simons & Olhede (2013)
%          4 SIMULOSL - single-field as in Simons & Olhede (2013)
% Hk       The Fourier-domain data, e.g. from SIMULOSL
% thhat    The evaluated scaled parameter vector
% scl      The scale factors
% params   The structure with the fixed parameters from the experiment
% ah       A triplet of axis handles [defaulted]
% pertur   A flag identifying "bad" results for axis scaling [default 0]
%          The following options only when calling from OLHEDESIMONS5
% th0      The true vector
% covX     The parameter covariance that you want quoted
% E,v      Young's modulus and Poisson's ratio, should you require them
%
% OUTPUT:
%
% cb,xl    Axis handles
% t,tt     Title and subtitle handles
%
% NOTE: 
%
% Watch what happens when you've filtered and check the warnings related
% to blurring, filtering, etc.
%
% SEE ALSO:
%
% MLECHIPSDOSL, MLEOSL, EGGERS6
%
% Last modified by olwalbert-at-princeton.edu, 06/01/2025
% Last modified by fjsimons-at-alum.mit.edu, 06/01/2025

defval('pertur',0)
defval('th0',[])
defval('covX',[])
defval('E',[])
defval('v',[])

% Color bar orientation
cborien='hor';
cborien='vert';

% So now you have found a solution and you evaluate it at thhat
switch witsj
  case 1
   [~,~,~,~,k,~,~,Sb,Lb]=simulos(thhat.*scl,params,0);
  case 2
   [~,~,~,~,k,~,~,Sb,Lb]=simulros(thhat.*scl,params,0);
  case 3
   [~,~,~,~,k,~,~,Sb,Lb]=simulros0(thhat.*scl,params,0);
  case 4
   [~,~,~,k,~,Sb,Lb]=simulosl(thhat.*scl,params,0);
end

switch witsj
  case {1,2,3}
   % Invert this L (or Lbar, as the case may be) matrix
   % Note that L is lower-triangular!
   detLb=Lb(:,1).*Lb(:,3);
   invLb=[Lb(:,3) -Lb(:,2) Lb(:,1)]./repmat(detLb,1,3);
   detSb=[Sb(:,1).*Sb(:,3)-Sb(:,2).^2];
   invSb=[Sb(:,3) -Sb(:,2) Sb(:,1)]./repmat(detSb,1,3);
   % Check the inverse of the symmetric matrices and the Cholesky
   invcheck(invSb,detSb,Sb,10,1)
   invcheck(invSb,detSb,Sb,10,2)
   cholcheck(Lb,Sb,6,1)
   cholcheck(Lb,Sb,6,2)
   
   % Actually, now you should note that LKROS, LKROS0, LKOS, and LKOSL
   % really do this already, this is nothing more than the distribution of
   % the coefficients that go into their sum to make the likelihood
   % itself. And I just put those bits in there - will be able to take
   % them out from here later!
   
   % Remember it's X being chi-2_4 divided by 2
   % Multiply to obtain a variable which should follow the rules
   Zk=[invLb(:,1).*Hk(:,1) [invLb(:,2).*Hk(:,1)+invLb(:,3).*Hk(:,2)]];
   % And then form the product Zk^H*Zk but instead
   Xk0=hformos(1,Zk,[1 0 1]);
   % Same thing
   Xk=hformos(1,Hk,invSb);
   % This should be the degrees of freedom of the chi-squared of 2*Xk
   df=4;
  case 4
    % This only works if params.blurs is not Inf since then SIMULOSL changed
    if ~any(isnan(Lb))
        % Multiply to obtain a variable which should follow the rules 
        Zk=[Hk(:)./Lb(:)];
        Xk0=hformos(1,Zk,[1 0 1]);
        % Same thing
        Xk=abs(Hk(:)).^2./Sb(:);
        diferm(Xk,Xk0)
    else
        [Xk0,Xk]=deal(NaN);
    end
    % This should be the degrees of freedom of the chi-squared of 2*Xk
    df=2;
end

% And this should be the same thing again, except how it treats k=0
% Note that if you HAVE a solution already, you'll find the loglik you had
switch witsj
  case 1
   Xk1=-Lkos(k,thhat.*scl,params,Hk)-log(detSb);
   Lbar=loglios(thhat,params,Hk,k,scl);
  case 2
   Xk1=-Lkros(k,thhat.*scl,params,Hk)-log(detSb);
   Lbar=logliros(thhat,params,Hk,k,scl);
  case 3
   Xk1=-Lkros0(k,thhat.*scl,params,Hk)-log(detSb);
   Lbar=logliros0(thhat,params,Hk,k,scl);
  case 4
   [Lbar,~,~,momx,~,Lk]=logliosl(k,thhat.*scl,params,Hk,1);
    % This only works if params.blurs is not Inf since then SIMULOSL changed
    if ~any(isnan(Lb))
        Xk1=-Lk-log(Sb(~~k));
        % Check we're doing the same thing to tolerance, depending on whether
        % some prior codes put a NaN at zero wavenumber or got rid of the
        % zero-wavenumber values altogether; these checks will be removed
        difer(Xk(~isnan(Xk0))-Xk0(~isnan(Xk0)),9,[],NaN)
        difer(Xk(~~k)-Xk1,9,[],NaN)
    else
        % You cannot recompute it from the output of LOGLIOSL, you must redo it
        % It might be a reason to save LKOSL after all
        params.blurs=-1;
        % [S,kk]=maternosp(th0,params,1);
        % Same thing
        S=blurosy(th0,params,1);
        % Remember Sb that came out of SIMULOSL in this case was the periodogram
        % and not the expected periodogram, and it was for a new run not for the Hk
        % that were input here, hence we needed to recompute the expected periodogram
        Xk=abs(Hk).^2./S;
    end
   % The oldest way, using a since retired function LKOSL
   % Xkk1=-lkosl(k,thhat.*scl,params,Hk)-log(Sb);
   % difer(Xkk1(~~k)-Xk1,9,[],NaN)
end

switch witsj
 case {1,2,3}
  titst=sprintf('L =  %8.3f   ln(det[S]) = %8.3f',-Lbar,-mean(log(detSb(detSb>0))));
  varibal='Xo';
 case 4
  titst=sprintf('L =  %8.3f   ln(S) = %8.3f',-Lbar,-mean(log(Sb(Sb>0))));
  varibal='X';
end

% First make the wavenumbers, give the data size and the data length
% k=knums(params);

% Evaluate the likelihood
% disp(sprintf('The loglihood is %8.3f',Lbar))
% Craft some labels
xll=[0 3*2*df];
xlls=[xll(1):df:xll(2)];
lx=sqrt(length(Xk));
xlis=[0.5 lx/2 lx+0.5];
xstr2=sprintf('quadratic residual 2%s',varibal);
xstr=sprintf('quadratic residual %s',varibal);
cax=[0 3*df];

% If the graphics handle by this name exists
if ~exist('ah','var')
  ah=krijetem(subnum(1,3));
end

% Select non-ridiculous values
allg=~isinf(Xk);
% Should that perhaps be isnan, as in MLECHIPSDOSL

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(1))
% Take out the Inf which may occur at zero wavenumber
[bdens,c]=hist(2*Xk(allg),5*round(log(length(Xk(allg)))));
% Plot the histogram as a bar graph
bdens=bdens/indeks(diff(c),1)/length(Xk(allg));
bb=bar(c,bdens,1);
% Each of the below should be df/2
t(1)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
		   varibal,nanmean(Xk(allg)),...
		   varibal,nanvar(Xk(allg))));
set(bb,'FaceColor',grey)
hold on
% Plot the ideal chi-squared distribution
refs=linspace(0,max(2*Xk),100);
plot(refs,chi2pdf(refs,df),'Linew',1,'Color','k')
hold off

% Labeling and cosmetic adjustments
if pertur==0
  maxi=0.25;
  ylls=[0:0.1:maxi*(4/df)];
  ylim([0 maxi*(4/df)])
else
  ylls=[0 0.1 0.2 0.3];
  ylim([0 0.35*(4/df)])
end
xlim(xll)
xl(1)=xlabel(xstr2); 
yl(1)=ylabel('probability density');
axis square

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))
% Note that SOME people use a different parameterization (b vs 1/b)
% Note that gamma [df/2 2] is chi-squared [df]...
% Note that the try/catch provides the necessary version upgrade
try
  h=qqplot(2*Xk(allg),makedist('gamma','a',df/2,'b',2));
catch
  h=qqplot(2*Xk(allg),ProbDistUnivParam('gamma',[df/2 2])); 
end

axis image; box on
set(h(1),'MarkerEdge','k')  
set(h(3),'LineStyle','-','Color',grey)
% Extend the line to the full axis
hold on
xh=get(h(3),'xdata');
yh=get(h(3),'ydata');
h(4)=plot([xh(2) xll(2)],...
	  [yh(2) yh(2)+[yh(2)-yh(1)]/[xh(2)-xh(1)]*[xll(2)-xh(2)]]);
set(h(4),'LineS','-','Color',grey)
hold off
top(h(3),ah(2))
delete(get(ah(2),'ylabel'));
delete(get(ah(2),'title'));

% Take a look a the distribution of the residual moments
% See RB X, p. 51 about the skewness of a chi-squared - just sayin'.
% We don't change the number of degrees of freedom! If you have used
% twice the number, and given half correlated variables, you do not
% change the variance, that is the whole point. Unlike in FISHIOSL
% where you make an analytical prediction that does depend on the
% number and which therefore you need to adjust.

% Test for departure of chi-squaredness via the "magic" parameter which
% is of course, sort of, a sample variance, and normally distributed by
% the law of large numbers (with the postulated population mean
% subsituted). At any rate, magx should be close to nanvar(2X) above.
% This is the same as what comes out of LOGLIOSL etc
magx=nanmean([Xk(allg)-df/2].^2);
neem='\xi'; neem='s_X^2';
% Do the test whether you accept this as a good fit, or not
vr=8/length(k(~~k));
% Then use NORMTEST to ascertain the veracity... don't bother with the
% Nyquist wavenumbers, there will be very few, but take out the zero
[a,b,c]=normtest(magx,1,vr,0.05);

% disp(sprintf('NORMTEST %i %5.3f %i',a,b,round(c)))
if a==0; stp='accept'; else stp='reject'; end
t(2)=title(sprintf('%s =  %5.3f   8/K = %5.3f   %s   p = %5.2f',...
                   neem,magx,vr,stp,b));
% This we need for when we have multiple experiments
% t(2)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
%	           neem,nanmean(magx),...
%	           neem,nanvar(magx));
delete(get(ah(2),'xlabel'));
xlim(xll); ylim(xll)
xl(2)=xlabel(sprintf('predicted 2%s',varibal));
yl(2)=ylabel(sprintf('observed 2%s',varibal));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(3))
imagefnan([xlis(1) xlis(end)],[xlis(end) xlis(1)],...
	  reshape(Xk,params.NyNx),...
	  'gray',cax,[],1); axis image
% Whatever you do, take out the zero wavenumber, where you get a special point
t(3)=title(titst);
set(ah(3),'xtick',xlis,'XtickLabel',[-lx/2 0 lx/2],...
	  'ytick',xlis,'YtickLabel',[-lx/2 0 lx/2])
[cb,xcb]=addcb(cborien,cax,cax,'gray',df,1);
set(xcb,'string',xstr)

% Cosmetics
fig2print(gcf,'landscape')
set(ah(1),'xtick',xlls,'ytick',ylls)
set(ah(2),'xtick',xlls,'ytick',xlls)
longticks([ah cb])

% Any off-putting motions would have been due to using underscores in titles
set(ah(3),'position',...
 	  [getpos(ah(3),1) getpos(ah(3),2) getpos(ah(2),3) ...
 	   getpos(ah(2),4)])

if strcmp(cborien,'hor')
  movev(ah(3),-0.09)      
  movev(cb,-0.08)
  shrink(cb,0.9,1.3)
  movev(xcb,-1.25)
elseif strcmp(cborien,'vert')
  axes(ah(3))
  xl(3)=xlabel('wavenumber index'); 
  yl(3)=ylabel('wavenumber index'); 
  moveh(yl(3),params.NyNx(1)/30)

  axes(cb)
  moveh(cb,.075)
  set(cb,'yaxisl','r')
  shrink(cb,1.3,1)
  set(cb,'position',...
	 [getpos(cb,1) getpos(ah(3),2) getpos(cb,3) getpos(ah(3),4)])
  shrink(cb,1,1.27)
  moveh(xcb,15)
  moveh([ah cb],-.025)
end

% Cosmetics for the labels
set([yl(:); xl(:); xcb(:)],'FontSize',11)

% Cosmetics for the title
axes(ah(1))
movev(t(1),0.02)
axes(ah(2))
movev(t(2),0.1)
axes(ah(3))
delete(t(3))
delete(xcb)
t(3)=title(xstr);
movev(t(3),2.5)

% Stick the params here somewhere so we can continue to judge
movev([ah cb],-.1)
tt=ostitle(ah,params,[],length(thhat(:,1)));
movev(tt,.35)
% E.g. quote  the TRUTHS and the THEORETICAL standard deviation with
% which it can be known using the available data... as you wish... or not
% If run from MLEOSL5 you'll want to supply the proper values if you have them
[answ,answs]=osansw(thhat.*scl,covX,E,v);
ttt=supertit(ah,sprintf(answs,answ{:}));
movev(ttt,-4.25)

if ~isempty(E) && ~isempty(v)
    movev(tt,-0.25)
    % But ALSO show the distance between the TRUTH and the estimate for
    % the effective elastic thickness Te in km. Don't forget the transform
    % is nonlinear! Don't subtract before transforming.
    disp(sprintf('ABS Distance of estimate to truth is %5.3g km',...
		 (abs(DtoTe(thhat(1)*scl(1),E,v)-...
		      DtoTe(th0(1),E,v))/1000)))
end

% Output
varns={cb,xl,t,tt};
varargout=varns(1:nargout);
