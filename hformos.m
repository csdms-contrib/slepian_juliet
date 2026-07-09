function Xk=hformos(S,Hk,A,xver)
% Xk=HFORMOS(S,Hk,A,xver)
%
% Computes a certain quadratic form used by Olhede & Simons, namely
% the function (1/S) [Hk^H * A * Hk] where ^H is the Hermitian transpose
%
% INPUT:
%
% S     A spectral density, which will be inverted
% Hk    A complex 2-column matrix of Fourier-domain observations, [Nx2], or
%       a complex 1-column vector of Fourier-domain observations, [Nx1]
% A     A wavenumber-dependent SYMMETRIC matrix, e.g. Tinv, [Nx3], or
%       a wavenumber-dependent vector, [Nx1]
% xver  1 or 2 Extra verification for the univariate case, when A=1
%
% OUTPUT:
%
% Xk   A column vector with the wavenumbers unwrapped, [N]
%
% NOTE: 
%
% In our formalism, this should result in a real-valued function. 
%
%
% EXAMPLE:
%
% hformos('demo1') % compare the uncorrelated residuals (Xk) to a similar 
%                  % quadratic form that includes wavenumber interaction,
%                  % and assess the chi-squaredness of their distributions
% hformos('demo2') % generate an anisotropic field to study the distribution of
%                  % the residual term Xk
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2023
% Last modified by olwalbert-at-princeton.edu, 05/21/2025

defval('A',1)
defval('xver',1)

if ~isstr(S)
    % The Inf's at k=0 get turned into NaN's here
    Xk=    1./S.*[  conj(Hk(:,1)).*A(:,1).*Hk(:,1)];
    if size(A,2)==3 && size(Hk,2)==2
        Xk=Xk+1./S.*[2*real(conj(Hk(:,2)).*A(:,2).*Hk(:,1))+...
		     conj(Hk(:,2)).*A(:,3).*Hk(:,2)];
    end

    % The 2*real is a shortcut to write the sum of the cross terms which
    % should cancel any imaginary parts... it's a magnitude after all
    % Still may have a tiny imaginary part, get rid of it if it is small 
    Xk=realize(Xk);

    if xver==1 || xver==2
        if A==1
            diferm(Xk,abs(Hk).^2./S,7);
        end
    end
elseif strcmp(S,'demo1')
  % Create visualizations of a residual term that incorporates wavenumber 
  % interaction. Equations are in the subplot titles (roughly).
  % Set the grid parameters (less that 120 by 120) and the Matern model
  % parameters
  p=[]; p.NyNx=[48 48]; p.dydx=[1 1]; p.blurs=-1; p.taper=1;
  sp=p; sp.blurs=Inf;
  % Note that the behavior of the distribution of the residuals for the
  % uncorrelated and correlated case depends on theta, particularly nu (best
  % when nu==1)
  %th=[1 0.5 1];
  th=[2 1.0 3.5];

  % We will be making comparisons between the uncorrelated X(k) and a 
  % correlated version through the terms involved in their calculation: H(k),
  % \bar{S}(k), and <H(k),H(k')>, which is a term from the Isserlis
  % calculations of COVGAMMIOSL that requires C(y). 

  % Simulate some correlated (blurs=Inf) spectral data, H(k). Since we need to 
  % be careful with tracking the organization of k and k', we will ask for the 
  % wavenumbers, too.
  [~,~,~,k,Hk]=simulosl(th,sp);
  HkHkp=Hk*Hk';

  % Create X(k) and \bar{S}(k) for comparison later
  Sb=blurosy(th,p,0,'efd',[0 0]);
  Sb=realize(fftshift(v2s(Sb,p))); Sb=Sb(:);
  Xk=hformos(Sb,Hk,[],1);

  % Syntax copied over from my covgammiosl method. I invert <H(k),H(k')>
  % (EJJht_out) and use a similar instrumentation of organizing the [ MN x MN ] 
  % matrix in the third dimension as [ M x N x MN ] so that I can multiply the 
  % rows of the [ M x N ] H(k) with the second dimension of the reshaped 
  % <H(k),H(k')>, which represent all ks and a specific k' in each of the MN pages

  %%% BEGIN COPY FROM COVGAMMIOSL, METHOD 2 (DFTMTX)
  dydx=p.dydx;NyNx=p.NyNx;                                        
  ys=0:NyNx(1)-1;xs=0:NyNx(2)-1;                                            
  [grdx,grdy]=meshgrid(xs(:),ys(:));                                        
  grdx=grdx(:).*dydx(2);                                                    
  grdy=grdy(:).*dydx(1);                                                    
  lagx=grdx-grdx';                                                          
  lagy=grdy-grdy';                                                          
  dxxp=sqrt(lagx.^2+lagy.^2);    
  Cmn=maternosy(dxxp,th);
  if numel(p.taper)>1
    Cmn=Cmn.*Tx(:).*Tx(:)';
  end
  Cm_n_mn=reshape(Cmn,[NyNx(1) NyNx(2) prod(NyNx)]);
  % The DFT matrix (mxm)
  Um=dftmtx(size(Cm_n_mn,1));
  % The DFT matrix (nxn)
  Un=dftmtx(size(Cm_n_mn,2));
  % Form the shared inner term of both EJJht and EJJt
  % This is like taking fft2(Cm_n_mn(:,:,i),[],2) for the ith slice
  EJJht_in=reshape(permute(...
            tensorprod(Um,tensorprod(Un,Cm_n_mn,1,2),1,2),...
           [3 1 2]),NyNx(1),NyNx(2),prod(NyNx));
           EJJht_out=conj(tensorprod(Um,tensorprod(Un,conj(EJJht_in),1,2),1,2));
  EJJht_out=reshape(permute(EJJht_out,[3 1 2]),[prod(NyNx) prod(NyNx)]);
  %%% END COPY FROM COVGAMMIOSL, METHOD 2 (DFTMTX)

  % We know that the diagonal of |EJJht_out|^2 is an unwrapping of the \bar{S}(k)
  % that we calculate from BLUROSY using blurs=-1 (exact blurring, uncorrelated)
  % from the study in COVAGMMIOSL DEMO4

  % Scale, calculate the inverse, and reshape following the same pattern as 
  % before
  expHkHkp=EJJht_out./(4*pi^2*prod(NyNx));
  em_n_mn=reshape(inv(expHkHkp),[NyNx(1) NyNx(2) prod(NyNx)]);
  % Be very careful that the wavenumbers are ordered correctly. If we end up
  % using this for calculations later, it will be worth confirming this with
  % tests many times over for rectangular grids with (un)symmetric tapers

  % Organize the data vector on the grid with zero wavenumber in the (1,1) spot
  mHk=fftshift(v2s(Hk,p));
  % This is the same multiplication pattern as above from COVGAMMIOSL, 
  % which we know works correctly
  rxkkp=permute(tensorprod(mHk',tensorprod(mHk,em_n_mn,1,2),1,2),[3 1 2]);
  Xkkp=(reshape((rxkkp),[prod(NyNx) prod(NyNx)]));

  % Would we ever want to average over k (mean(rxkkp,1)) or k' (mean(rxkkp,2))?
  % % mrxkkp=squeeze(mean(rxkkp,1));
  %
  % Would we ever want to reorganize k and k' to look in at the residual from a
  % different perspective? No option seems to offer much visual intuition
  % % xkkp=fftshift(reshape(fftshift(rxkkp),[prod(NyNx) prod(NyNx)]));
  
  % Grab the diagonal of <H(k),H(k')> and reshape it [ M x N ]. This is what I
  % compare to \bar{S}(k). 
  dexpHkHkp=fftshift(v2s(diag(realize(expHkHkp)),p));
  % Grab the diagonal of H(k)'[<H(k),H(k')>]^{-1}H(k') and reshape to [M x N]
  % for comparison with X(k)
  % dXkkp=realize(v2s(diag(Xkkp),p));

  % During the 5/20 meeting I was showing this for the diagonal
  dXkkp=realize(v2s(diag(HkHkp),p)./dexpHkHkp); 

  % The magnitude and phase of H(k)'[<H(k),H(k')>]^{-1}H(k')
  aXkkp=abs(Xkkp);
  pXkkp=angle(Xkkp);

  % Now for the figures
  fs=13; % font-size for the latex interpreted subplot titles
  
  % Figure 1: visualize the k==k' terms (i.e., |H(k)|^2, \bar{S}(k), X(k), or
  % the main diagonal of their correlated counterparts that I calculate above)
  figure(1)
  subplot(331); imagesc(v2s(Xk,p)); cax=caxis;
  xlabel('k_x'); ylabel('k_y');
  title('$X(k)=\frac{|H(k)|^2}{\bar{S}}$','Interpreter','latex','FontSize',fs)
  subplot(332); imagesc(dXkkp); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  %title('$\mathrm{diag}(X(k,k''))=\mathrm{diag}(\frac{H(k)H(k'')}{<H(k),H(k'')>})$',...
  title('$\mathrm{diag}(X(k,k''))=\mathrm{diag}(H(k)'' [<H(k),H(k'')>]^{-1} H(k''))$',...
        'Interpreter','latex','FontSize',fs)
  subplot(333); imagesc(v2s(Xk,p)./dXkkp); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$\frac{Xk}{\mathrm{diag}(X(k,k''))}$','Interpreter','latex','FontSize',fs)

  subplot(334); imagesc(v2s(Sb,p)); cax=caxis;
  xlabel('k_x'); ylabel('k_y');
  title('$\bar{S}(k)$','Interpreter','latex','FontSize',fs)
  subplot(335); imagesc(dexpHkHkp); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$\mathrm{diag}(<H(k),H(k'')>)$','Interpreter','latex','FontSize',fs)
  subplot(336); imagesc(v2s(Sb,p)./dexpHkHkp); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$\frac{\bar{S}}{\mathrm{diag}(<H(k),H(k'')>)}$',...
        'Interpreter','latex','FontSize',fs)

  subplot(337); imagesc(abs(v2s(Hk,p)).^2); cax=caxis;
  xlabel('k_x'); ylabel('k_y');
  title('$|H(k)|^2$','Interpreter','latex','FontSize',fs)
  subplot(338); imagesc(abs(v2s(diag(HkHkp),p))); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$|\mathrm{diag}(H(k)''H(k''))|$','Interpreter','latex','FontSize',fs)
  subplot(339); imagesc(abs(v2s(Hk,p)).^2./abs(v2s(diag(HkHkp),p))); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$\frac{|H(k)|^2}{|\mathrm{diag}(H(k)H(k''))|}$',...
        'Interpreter','latex','FontSize',fs)

  % Figure 2: visualize the magnitude and phase of the correlated terms whose
  % diagonals we just inspected in Figure 1
  figure(2)
  a1(1)=subplot(321); imagesc(aXkkp);
  xlabel('k_x'); ylabel('k_y');
  %title('$|X(k,k'')|=\frac{|HkHk''|}{|<Hk,Hk''>|}$',...
  title('$|(X(k,k'')|=|H(k)'' [<H(k),H(k'')>]^{-1} H(k'')|$',...
        'Interpreter','latex','FontSize',fs)
  a2(1)=subplot(322); imagesc(pXkkp);
  xlabel('k_x'); ylabel('k_y');
  %title('$\mathrm{phase}(X(k,k''))=\frac{\mathrm{phase}(HkHk'')}{\mathrm{phase}(<Hk,Hk''>)})$',...
  title('$\mathrm{phase}(X(k,k''))=\mathrm{phase}(H(k)'' [<H(k),H(k'')>]^{-1} H(k''))$',...
        'Interpreter','latex','FontSize',fs)

  a1(2)=subplot(323); imagesc(abs(expHkHkp));
  xlabel('k_x'); ylabel('k_y');
  title('$|<H(k),H(k'')>|$','Interpreter','latex','FontSize',fs)
  a2(2)=subplot(324); imagesc(angle(expHkHkp)); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$\mathrm{phase}(<H(k),H(k'')>)$','Interpreter','latex','FontSize',fs)

  a1(3)=subplot(325); imagesc(abs(HkHkp)); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$|H(k)''H(k'')|$','Interpreter','latex','FontSize',fs)
  a2(3)=subplot(326); imagesc(angle(HkHkp)); %caxis(cax)
  xlabel('k_x'); ylabel('k_y');
  title('$\mathrm{phase}(H(k)''H(k''))$','Interpreter','latex','FontSize',fs)

  % Phase is always too much to look at, tone it down with grayscale
  for ind=1:3
    colormap(a1(ind),'parula')
    colormap(a2(ind),'gray')
  end
  
  % Now we really should take a look at the distribution of X(k), diag(X(k,k'))
  % -- meaning the diagonal of the correlated counterpart to X(k), 
  % ``X(k,k')'' = H(k)'[<H(k),H(k')>]^{-1}H(k') , as well as
  % its magnitude and phase. While we do not necessarily want to express the
  % distribution of the correlated residuals in terms of chi-squared, I plot 
  % reference curves for 1, 2, and 3 degrees of freedom anyway, and incorporate 
  % the analysis of residuals that we do understand for uncorrelated residuals
  % (|H(k)|^2/\bar{S} ~ \chi^2_2 / 2 . The distribution of the phase seems to
  % always be uniform.

  % Subplots below are copied from MLECHI{PLOS,PSDOSL} with minimal modifications
  
  % Bounds of X residuals to show (functions as color axis in panel 3)
  df=2;
  boundsX0=[0 3*df];
  % Bounds of 2X residuals to how (functions as axis bounds in panel 1)
  bounds2X=2.*boundsX0;
  % The degrees of freedom that we will make reference plots for and some colors
  dfs=1:1:3;
  cs=parula(numel(dfs)+2);

  figure(3)
  ah=krijetem(subnum(2,2));
  % Xk first
  axes(ah(1))
  varibal='X(k)'; xstr2v=sprintf('quadratic residual 2%s',varibal);
  allg=~isinf(Xk);
  % Bins to hit the grid lines
  binWidth=df/4;
  binos=(binWidth/2)*[1:2:bounds2X(2)*2/binWidth+1];
  [bdens,c]=hist(2*Xk(allg),binos);
  % [bdens,c]=hist(2*Xk(allg),5*round(log(length(Xk(allg)))));
  % Plot the histogram as a bar graph
  bdens=bdens/indeks(diff(c),1)/length(Xk(allg));
  bb=bar(c,bdens,1,'HandleVisibility','off');
  set(bb,'FaceC',grey)
  hold on
  % Plot some reference chi-squared distributions with different dof
  refs=linspace(0,max(bounds2X),100);
  % refs=linspace(0,max(2*Xk),1000);
  for ind=1:numel(dfs)
    p1(ind)=plot(refs,chi2pdf(refs,dfs(ind)),'Linew',1,'Color',cs(ind,:),...
         'DisplayName',sprintf('dof=%i',dfs(ind)));
  end
  set(p1(2),'Color','k','LineWidth',2)
  lg(1)=legend();set(lg(1),'Location','north')
  hold off
  % Labeling and cosmetic adjustments
  xlim(bounds2X)
  xl1(1)=xlabel(xstr2v);
  ylim([0 max(bdens)*1.05])
  % ylim([0 0.35*(4/df)])
  yl1(1)=ylabel('probability density');
  % Each of the below should be df/2
  t(1)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
             varibal,nanmean(Xk(allg)),...
             varibal,nanvar(Xk(allg))));

  % diag(X(k,k')) second
  axes(ah(2))
  varibal='diag(X(k,k''))'; xstr2v=sprintf('quadratic residual 2%s',varibal);
  dXkkp=dXkkp(:);
  allg=~isinf(dXkkp);
  % Bins to hit the grid lines
  binWidth=df/4;
  binos=(binWidth/2)*[1:2:bounds2X(2)*2/binWidth+1];
  [bdens,c]=hist(2*dXkkp(allg),binos);
  % [bdens,c]=hist(2*dXkkp(allg),5*round(log(length(dXkkp(allg)))));
  % Plot the histogram as a bar graph
  bdens=bdens/indeks(diff(c),1)/length(dXkkp(allg));
  bb=bar(c,bdens,1,'HandleVisibility','off');
  set(bb,'FaceC',grey)
  hold on
  % Plot some reference chi-squared distributions with different dof
  refs=linspace(0,max(bounds2X),100);
  % refs=linspace(0,max(2*dXkkp),1000);
  dfs=1:1:3;
  cs=parula(numel(dfs)+2);
  for ind=1:numel(dfs)
    p2(ind)=plot(refs,chi2pdf(refs,dfs(ind)),'Linew',1,'Color',cs(ind,:),...
         'DisplayName',sprintf('dof=%i',dfs(ind)));
  end
  set(p2(2),'Color','k','LineWidth',2)
  lg(2)=legend();set(lg(2),'Location','north')
  hold off
  % Labeling and cosmetic adjustments
  xlim(bounds2X)
  xl1(2)=xlabel(xstr2v);
  ylim([0 max(bdens)*1.05])
  yl1(2)=ylabel('probability density');
  % Each of the below should be df/2
  t(2)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
             varibal,nanmean(dXkkp(allg)),...
             varibal,nanvar(dXkkp(allg))));

  % The magnitude of X(k,k') third
  axes(ah(3))
  varibal='|X(k,k'')|'; xstr2v=sprintf('quadratic residual 2%s',varibal);
  aXkkp=aXkkp(:);
  allg=~isinf(aXkkp);
  % Treat it just like Xk, with the exception that now we have complex values
  % Bins to hit the grid lines
  binWidth=df/4;
  binos=(binWidth/2)*[1:2:bounds2X(2)*2/binWidth+1];
  [bdens,c]=hist(2*aXkkp(allg),binos);
  % [bdens,c]=hist(2*aXkkp(allg),5*round(log(numel(aXkkp(allg)))));
  % Plot the histogram as a bar graph
  bdens=bdens/indeks(diff(c),1)/numel(aXkkp(allg));
  bb=bar(c,bdens,1,'HandleVisibility','off');
  set(bb,'FaceC',grey)
  hold on
  % Plot some reference chi-squared distributions with different dof
  refs=linspace(0,max(bounds2X),100);
  % refs=linspace(0,max(2*aXkkp),1000);
  dfs=1:1:3;
  cs=parula(numel(dfs)+2);
  for ind=1:numel(dfs)
    p3(ind)=plot(refs,chi2pdf(refs,dfs(ind)),'Linew',1,'Color',cs(ind,:),...
         'Linew',1,'Color',cs(ind,:),...
         'DisplayName',sprintf('dof=%i',dfs(ind)));
  end
  set(p3(2),'Color','k','LineWidth',2)
  lg(3)=legend();set(lg(3),'Location','east')
  hold off
  % Labeling and cosmetic adjustments
  xlim(bounds2X)
  xl1(3)=xlabel(xstr2v);
  ylim([0 max(bdens)*1.05])
  yl1(3)=ylabel('probability density');
  t(3)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
             varibal,nanmean(aXkkp(allg)),...
             varibal,nanvar(aXkkp(allg))));

  % The phase of X(k,k') last
  axes(ah(4))
  varibal='phase(X(k,k''))';
  pXkkp=pXkkp(:);
  allg=~isinf(pXkkp);
  [bdens,c]=hist(2*pXkkp(allg),5*round(log(numel(pXkkp(allg)))));
  % Plot the histogram as a bar graph
  bdens=bdens/indeks(diff(c),1)/numel(pXkkp(allg));
  bb=bar(c,bdens,1,'HandleVisibility','off');
  set(bb,'FaceC',grey)
  xl1(4)=xlabel('radians');
  ylim([0 max(bdens)*2])
  yl1(4)=ylabel('probability density');
  %xlim([0 prctile(pXkkp(allg),99.7)])
  hold off
  t(4)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
             varibal,nanmean(pXkkp(allg)),...
             varibal,nanvar(pXkkp(allg))));


  % Prepare an overlay axis for the quantile-quantile plots for subplot 1--3
  resterm={Xk,dXkkp,aXkkp,pXkkp};
  for ind=1:3
    allg=~isinf(resterm{ind});
    ah2(ind)=laxis(ah(ind),0,0);
    axes(ah2(ind))
    % Obtain qq-plot data
    % Note that the try/catch provides the necessary version upgrade
    try
      h=qqplot(2*resterm{ind}(allg),makedist('gamma','a',df/2,'b',2));
    catch
      h=qqplot(2*resterm{ind},ProbDistUnivParam('gamma',[df/2 2]));
    end
    hx=get(h,'Xdata'); hx=hx{1};
    hy=get(h,'ydata'); hy=hy{1};
    delete(h)
    % Make the qq-plot
    qq0=plot(bounds2X,bounds2X,'k'); hold on
    qq=plot(hx,hy,'LineS','none','Marker','o','MarkerF','r',...
    	'MarkerE','r','MarkerS',2);
    % More cosmetic adjustments
    set(ah(ind),'box','off')
    set(ah2(ind),'xlim',get(ah(ind),'xlim'),'yaxisl','r','box','off',...
    	   'xaxisl','t','color','none','ylim',get(ah(ind),'xlim'))
    % Add labels and tickmarks
    ylr(ind)=ylabel('quantile-quantile prediction');
    tkmks=bounds2X(1):2*df:bounds2X(2);
    set([ah(ind) ah2(ind)],'xtick',tkmks)
    set(ah2(ind),'XTickLabel',repmat({''},numel(tkmks)))
    set(ah2(ind),'ytick',tkmks)
    % Test for departure of chi-squaredness... copied from mlechiplos
    magx=nanmean([resterm{ind}(allg)-df/2].^2);
    neem='\xi'; neem='s_X^2';
    % Do the test whether you accept this as a good fit, or not
    vr=8/length(k(~~k));
    [a,b,c]=normtest(magx,1,vr,0.05);
    % disp(sprintf('NORMTEST %i %5.3f %i',a,b,round(c)))
    if a==0; stp='accept'; else stp='reject'; end
    t2(ind)=title(sprintf('%s =  %5.3f   8/K = %5.3f   %s   p = %5.2f',...
                  neem,magx,vr,stp,b));
    movev(t2(ind),0.3)
  end

   % % A SAMPLING APPROACH THAT SHOULD BE COMPARABLE TO THE MAIN CALCULATIONS
   % % ABOVE
   % 
   % [Hx,~,~,k,Hk]=simulosl(th,sp);
   % HkHkp=Hk*Hk';
   % Sb=blurosy(th,p);
   % Xk=hformos(Sb,Hk,[],1);
   % numreals=5000;
   % Hks=zeros(prod(p.NyNx),numreals);
   % ks=zeros(prod(p.NyNx),numreals);
   % parfor ind=1:numreals
   %   [Hx,~,~,k,Hk]=simulosl(th,sp);
   %   Hks(:,ind)=Hk;
   %   ks(:,ind)=k(:);
   % end
   %
   % % For testing the unwrapping:
   % % kkp=k(:)*k(:)';
   % % expkkp0=cov(ks');
   % % Compare kkp and expkkp
   % % diferm(kkp-expkkp0)
   % % expkkp1=ks*ks'./numreals;
   % % diferm(kkp-expkkp1)
   % % clf;subplot(121);imagesc(kkp./expkkp0);subplot(122);imagesc(kkp./expkkp1)
   % 
   % % Based on how the wavenumbers unwrap, it is better to take the product of 
   % % the data vectors for their covariance instead of using cov so that we 
   % % can understand their orientation
   % % expHkHkp=cov(Hks');
   %
   % expHkHkp=Hks*Hks'./numreals;
   % % The ratio we are after, calculated not quite the same as mth 2
   % Xkkp=HkHkp./expHkHkp;
   % aXkkp=abs(Xkkp);
   % pXkkp=angle(Xkkp);
   % dXkkp=v2s(diag(aXkkp),p);
   % dexpHkHkp=v2s(diag(expHkHkp),p);
   % %

   % % Some additional calculation paths considered on 5/20
   %
   % Xkkp=zeros(p.NyNx);
   % for knd=1:numel(Hk)
   %   for pnd=1:numel(Hk)
   %   %  Xkkp(knd,pnd)=circshift(Hk,knd-1)'*inv(expHkHkp)*circshift(Hk,pnd-1);
   %   end
   % end
   %
   % % Frederik's Dahlen and Simons 2008 idea: f'Ff=trace(ff'F)
   % em_n_mn=reshape(inv(expHkHkp),[NyNx(1) NyNx(2) prod(NyNx)]);
   % Xkkp=reshape(tensorprod(v2s(Hk,p)'*v2s(Hk,p),em_n_mn,1,2),[prod(NyNx) prod(NyNx)]);
   % trace(Xkkp)
   %
   % Frederik's note from meeting with Arthur 5/21 and previous conversations: 
   % blurred, correlated likelihood function?
   %
   % Sum or average over k or k'?
   % %
elseif strcmp(S,'demo2')
   % What structure do the residuals of an anisotropic Matern field exhibit when
   % modeled as an isotropic Matern process?
   
   % Simulate an anisotropic stationary Matern field by manipulating the dydx,
   % to be replaces with a call to anisimulosl in the future
   th=[2 1.2 5]; p.NyNx=[128 128]; p.dydx=[1 4]; p.blurs=Inf;
   [Hx,~,~,k,Hk]=simulosl(th,p);

   % Estimate using the isotropic framework
   p.blurs=-1; p.dydx=[1 1];
   [thhat,~,~,scl]=mleosl(Hx,[],p);
   Sb=blurosy(thhat.*scl,p);
   Xk=hformos(Sb,Hk);

   % Display the field, residuals, and distribution of the residuals with a QQ
   % plot and the chi-squaredness test
   fs=12;
   figure();
   subplot(321)
   imagesc(v2s(Hx,p))
   xlabel('x1'); ylabel('x2')
   longticks
   title('H(x)','Interpreter','latex','FontSize',fs)

   subplot(322)
   imagesc(abs(v2s(Hk,p)).^2)
   xlabel('k1'); ylabel('k2')
   longticks
   title('$|H(k)|^2$','Interpreter','latex','FontSize',fs)

   subplot(323)
   imagesc(v2s(Sb,p))
   xlabel('k1'); ylabel('k2')
   longticks
   title('$\bar{S}(k)$','Interpreter','latex','FontSize',fs)

   subplot(324)
   imagesc(v2s(Xk,p))
   xlabel('k1'); ylabel('k2')
   longticks
   title('$X(k)=|H(k)|^2/\bar{S}(k)$','Interpreter','latex','FontSize',fs)

   ah(5)=subplot(313);
   df=2;
   boundsX0=[0 3*df];
   % Bounds of 2X residuals to how (functions as axis bounds in panel 1)
   bounds2X=2.*boundsX0;
   varibal='X(k)'; xstr2v=sprintf('quadratic residual 2%s',varibal);
   allg=~isinf(Xk);
   % The degrees of freedom that we will make reference plots for and some colors
   dfs=1:1:3;
   cs=parula(numel(dfs)+2);
   % Bins to hit the grid lines
   binWidth=df/4;
   binos=(binWidth/2)*[1:2:bounds2X(2)*2/binWidth+1];
   [bdens,c]=hist(2*Xk(allg),binos);
   % Plot the histogram as a bar graph
   bdens=bdens/indeks(diff(c),1)/length(Xk(allg));
   bb=bar(c,bdens,1,'HandleVisibility','off');
   set(bb,'FaceC',grey)
   hold on
   % Plot some reference chi-squared distributions with different dof
   refs=linspace(0,max(bounds2X),100);
   for ind=1:numel(dfs)
     p1(ind)=plot(refs,chi2pdf(refs,dfs(ind)),'Linew',1,'Color',cs(ind,:),...
          'DisplayName',sprintf('dof=%i',dfs(ind)));
   end
   set(p1(2),'Color','k','LineWidth',2)
   lg(1)=legend();set(lg(1),'Location','north')
   hold off
   % Labeling and cosmetic adjustments
   xlim(bounds2X)
   xl1(1)=xlabel(xstr2v);
   ylim([0 max(bdens)*1.05])
   % ylim([0 0.35*(4/df)])
   yl1(1)=ylabel('probability density');
   % Each of the below should be df/2
   t(1)=title(sprintf('m(%s) =  %5.3f   v(%s) =  %5.3f',...
              varibal,nanmean(Xk(allg)),...
              varibal,nanvar(Xk(allg))));
   % The QQ plot 
   allg=~isinf(Xk);
   ah2=laxis(ah(5),0,0);
   axes(ah2)
   % Obtain qq-plot data
   % Note that the try/catch provides the necessary version upgrade
   try
     h=qqplot(2*Xk(allg),makedist('gamma','a',df/2,'b',2));
   catch
     h=qqplot(2*Xk,ProbDistUnivParam('gamma',[df/2 2]));
   end
   hx=get(h,'Xdata'); hx=hx{1};
   hy=get(h,'ydata'); hy=hy{1};
   delete(h)
   % Make the qq-plot
   qq0=plot(bounds2X,bounds2X,'k'); hold on
   qq=plot(hx,hy,'LineS','none','Marker','o','MarkerF','r',...
   	'MarkerE','r','MarkerS',2);
   % More cosmetic adjustments
   set(ah(5),'box','off')
   set(ah2,'xlim',get(ah(5),'xlim'),'yaxisl','r','box','off',...
    	   'xaxisl','t','color','none','ylim',get(ah(5),'xlim'))
    % Add labels and tickmarks
    ylr=ylabel('quantile-quantile prediction');
    tkmks=bounds2X(1):2*df:bounds2X(2);
    set([ah(5) ah2],'xtick',tkmks)
    set(ah2,'XTickLabel',repmat({''},numel(tkmks)))
    set(ah2,'ytick',tkmks)
    % Test for departure of chi-squaredness... copied from mlechiplos
    magx=nanmean([Xk(allg)-df/2].^2);
    neem='\xi'; neem='s_X^2';
    % Do the test whether you accept this as a good fit, or not
    vr=8/length(k(~~k));
    [a,b,c]=normtest(magx,1,vr,0.05);
    % disp(sprintf('NORMTEST %i %5.3f %i',a,b,round(c)))
    if a==0; stp='accept'; else stp='reject'; end
    t2=title(sprintf('%s =  %5.3f   8/K = %5.3f   %s   p = %5.2f',...
                  neem,magx,vr,stp,b));
    movev(t2,0.3)
end

