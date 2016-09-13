function varargout=simulosl(th0,params,xver,varargin)
% [Hx,th0,params,k,Hk,Sb,Lb,gane,miy]=SIMULOSL(th0,params,xver)
% [gane,miy,ax,legs]=SIMULOSL('demo4',th0,params,xver,ptype,N,rindj,Nspm,ah,gane,miy)
%
% Simulates single-field data in the Matern form. 
%
% INPUT:
%
% th0      The true parameter vector with elements:
%          th0(1)=s2   The first Matern parameter, aka sigma^2 
%          th0(2)=nu   The second Matern parameter 
%          th0(3)=rho  The third Matern parameter 
% params   A structure with constants that are (assumed to be) known:
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact BLUROSY procedure
%          kiso  wavenumber beyond which we are not considering the spectrum
%          quart 1 quadruple, then QUARTER
%                0 size as is, watch for periodic correlation behavior
% xver     1 for extra verification, 0 if not needed
% Only for 'demo4', which is used by EGGERS4
% ptype    'mle' or 'poor' for 'demo4' only [default: 'poor']
% N         Maximum size that we will be trying [default: 128]
% rindj     Steps of sizes that are being tried... [default: 2:2:N]
% Nspm      This many experiments to each of the processors [default: 5]
% ah        Axis handles so we can see if it's a multipanel figure or not
%           and also, the right panel gets the ylims from the left panel
% gane     Numbers for the axis equalization procedure in 'demo4'
% miy      Numbers for the axis equalization procedure in 'demo4'
%
% OUTPUT:
%
% Hx       Real matrix of spatial-domain observations [m], see Hk
% th0      The true parameter vector pertaining to this data
% params   The structure with the known knowns, see above
% k        Wavenumber(s) suitable for the data sets returned [rad/m]
% Hk       A complex matrix of Fourier-domain observations, namely
%          final surface topography [m]
% Sb       The spectral matrix that you've used in this process
% Lb       The Cholesky decomposition of the spectral matrix which you
%          might use to evaluate the fit later on, in which case you
%          don't need any of the previous output
% Only for 'demo4':
% gane     Numbers for the axis equalization procedure
% miy      Numbers for the axis equalization procedure
% ax       Handle to the extra axis
% legs     Handles to the objects that you're most likely to want legends for 
%
% EXAMPLE:
%
% [Hx,th0,p]=simulosl; p.blurs=-1; imagesc(reshape(simulosl(th0,p),p.NyNx))
%
% simulosl('demo1') % Just a little plot to illustrate the default behavior
% simulosl('demo2') % Only for symmetry with SIMULOS, SIMULROS, etc 
% simulosl('demo3') % Plots a couple of likelihoods to test LKOSL also
% simulosl('demo4') % Naive and mle variance and their biases, used by EGGERS4
% simulosl('demo5',th0,p) % A really poor space-domain covariance estimator
% simulosl('demo6',th0,p) % A better space-domain covariance estimator
%
% SEE ALSO:
%
% MLEOSL, LOADING, SIMULOS, EGGERS1, EGGERS2
%
% Last modified by fjsimons-at-alum.mit.edu, 05/03/2016

% Here is the true parameter vector and the only variable that gets used 
defval('th0',[1e6 2.5 20e3]);

% If not a demo...
if ~isstr(th0)
  % Supply the needed parameters, keep the givens, extract to variables
  fields={               'dydx','NyNx','blurs','kiso','quart'};
  defstruct('params',fields,...
	    {                      [10 10]*1e3,[128 128]/2,-1,NaN,1});
  struct2var(params)

  % Here is the extra verification parameter
  defval('xver',1)
  if xver==1
    % Dump to screen
    osdisp(th0,params)
  end

  % If you're going to be quartering, must be quadrupling first!
  % Actually, THAT is the question - or do we address this using tospace
  % and tospec, using zero-padding, i.e. work with the original wavevecs?
  if quart==1
    params.NyNx=NyNx*2;
    struct2var(params)
  end

  % First make the wavenumbers, given the data size and the data length
  [k,dci,dcn]=knums(params);

  % This should make sense as the spacing in wavenumber domain
  dkxdky=2*pi./NyNx./dydx;

  % Now construct the whole-spectral matrix
  Z1=randgpn(k,dci,dcn);
  % This cramps the style. But still. Should I put in a zero at dci? 
  % disp(sprintf('Z1: mean %+6.3f ; stdev %6.3f',...
  %    mean(Z1(:)),std(Z1(:))))

  switch blurs
   case {0,1}
    % disp(sprintf('%s without blurring',upper(mfilename)))
    % Now make the spectral-spectral portion of the spectral matrix
    S=maternos(k,th0);
    % Roll in the sqrt of the factored portion
    Lb=sqrt(S);
    % The spectral matrix in case you care but you don't
    Sb=S;
   otherwise
    % If I stay on the same k-grid, I'm really not doing any convolution at
    % all as you can see quickly. So by "suitably discretizing" the
    % convolutional operator we mean performing it on a highly densified
    % grid of which the target grid may be a subset. Doing this on the
    % same grid would be the "inverse crime" of not changing the grid at
    % all. Run Fk for this case to see it then would be a delta function

    if blurs>1
      % Blurs IS the refinement parameter; make new wavenumber grid
      % Then make the spectral-spectral portion of the spectral matrix
      % Then do the blurring and subsampling to original grid
      Sb=bluros(maternos(knums(params,1),th0),params,xver);
    else
      % Here is the alternative EXACT way of doing it, which does away
      % with the approximate convolutional refinement procedure
      Sb=blurosy(th0,params,xver);
    end

    % And then we do the Cholesky decomposition of that, explicitly
    Lb=sqrt(Sb);

    % Should make sure that this is real! Why wouldn't it be?
    Lb=realize(Lb);
  end
  
  % Blurred or unblurred, go on

  % And put it all together, unwrapped over k and over x
  Hk=Lb.*Z1(:);

  % Zero out the corners if so asked
  if any(~isnan(params.kiso))
    disp('SIMULOSL: kiso comes into play!')
    Hk(k>params.kiso)=0;
  end

  % And go to the space domain - unitary scaled transform
  % Watch the 2pi in MLEOSL
  Hx=(2*pi)*tospace(Hk,params);
  if xver==1
    % Check Hermiticity before transformation, absolute tolerance
    hermcheck(reshape(Hk,NyNx))
    % Check unitarity of the transform; relative tolerance
    diferm(Hk-tospec(Hx,params)/(2*pi),[],9-round(log10(mean(abs(Hk)))));
  end
  
  % Big it down
  if quart==1
    Hx=quarter(reshape(Hx,NyNx));
    Hx=Hx(:);
    % Undo what was done above
    params.NyNx=NyNx/2;
    struct2var(params)
    [k,dci,dcn]=knums(params);
    % Now do not forget that Hk, Sb, Lb etc are still on the doubled grid
    % This may have implications later on, that we choose to ignore now,
    % an example is if some of this output were to be passed onto LKOSL
  end
  
  % Return the output if requested
  varns={Hx,   th0,params,k,Hk,   Sb,Lb};
  varargout=varns(1:nargout);
elseif strcmp(th0,'demo1')
  % Pulls out the demo name
  svnm=th0
  % Does the simulation straight from this very own filename, all defaults
  [Hx,   th0,p,k,Hk   ]=feval(mfilename);
  struct2var(p)

  clf
  kelicol
  [ah,ha]=krijetem(subnum(2,2));
  delete(ah(3:4))
  ah=ah(1:2);
  [tl(1),cb(1),xc(1),xa(1)]=plotit(ah(1),Hx/1000,size(k),...
				   'Matern surface','field (%s)','km');

  % Cosmetics
  movev(ah,-.25)
  they=linspace(1,NyNx(1),5);
  thex=linspace(1,NyNx(2),5);
  spunkm=(NyNx-1).*dydx/1000;
  set(ah,'ylim',they([1 end])+[-1 1]/2,...
	 'xlim',thex([1 end])+[-1 1]/2,...
	 'ytick',they,...
	 'xtick',thex,...
	 'ytickl',-spunkm(1)/2+(they-1)*dydx(1)/1000,...
	 'xtickl',-spunkm(2)/2+(thex-1)*dydx(2)/1000)
  longticks([ah cb])

  % Plot the parameters here
  axes(ah(2))
  axis([-1 1 -1 1])
  nolabels(ah(2)); noticks(ah(2)); box on; axis image
  
  xof=-0.75;
  tx(1)=text(xof, 0.00,sprintf('%s = %12.3g','\sigma^2',th0(end-2)));
  tx(2)=text(xof,-0.25,sprintf('%s = %12.3g','\nu',     th0(end-1)));
  tx(3)=text(xof,-0.50,sprintf('%s = %12.3g','\rho',    th0(end)));

  fig2print(gcf,'portrait')
  figdisp([],svnm,[],1,'pdf')
elseif strcmp(th0,'demo2')
  % Try to get close to the example we had in 2000
  NyNx=[100 100];
  s2=0.001;
  nu=0.5;
  rho=30000;
  z2=15000;
  % Perform the simulation... not complete
  [Hx,   th0,k,Hk,   params]=simulosl;
elseif strcmp(th0,'demo3')
  % Input some trial parameters
  % Cannot use QUARTERING without adjusting Hk
  params.quart=1;
  % With a negative the data size can grow a lot
  params.blurs=-1;
  params.NyNx=randi(200,[1 2]);
  % Do the simulations
  [Hx,th0,params,k,Hk,Sb,Lb]=feval(mfilename,[],params); 

  % Take a good look
  clf
  ah=krijetem(subnum(1,2));
  % Must recompute the wavenumbers if you are to use Hk, it's a hack
  if params.quart==1
    params.NyNx=params.NyNx*2;
    k=knums(params);
  end

  % Check out these two different things and think about them
  Lk1=Lkosl(k,th0,params,Hk);
  axes(ah(1))
  imagesc(decibel(reshape(Lk1,params.NyNx)))
  axis image
  title(sprintf('The likelihood computed from Hk using %s','LKOSL'))
  xlabel(sprintf('Average as LOGLIOSL would have it %8.5g',-nanmean(Lk1)))

  % Must reset if you are going to provide your own computation
  if params.quart==1
    params.NyNx=params.NyNx/2;
    k=knums(params);
  end
  
  % Check out the difference (or ratio) between Hk and Hx and between Lk1
  % and Lk2 which should be nothing - and of course it is, watching (2pi)
  axes(ah(2))
  Hk2=tospec(Hx,params)/(2*pi);
  Lk2=Lkosl(k,th0,params,Hk2); 
  imagesc(decibel(reshape(Lk2,params.NyNx)))
  axis image
  title(sprintf('The likelihood computed from Hx using %s','LKOSL'))
  xlabel(sprintf('Average as LOGLIOSL would have it %8.5g',-nanmean(Lk2)))

  % Mark the size
  set(ah,'xtick',[1 params.NyNx(2)],'ytick',[1 params.NyNx(1)])
  set(ah,'clim',[0 3])
elseif strcmp(th0,'demo4')
  % This was written especially for EGGERS4!
  % Shift the inputs so you can default them but also supply them
  % No further inputs needed, but if you have them, you keep them
  % No inputs needed, but if you had them, you should use them in this order
  defval('params',[]); th0=params; clear params
  % Remember that for odd sample sizes 'rindj' 'blurs' better be odd also
  defval('xver',  []);   p  =xver; clear xver
  % More input if you so desire
  if nargin>3; ptype=varargin{1}; end; defval('ptype','poor')
  % Maximum size that we will be trying
  if nargin>4; N=varargin{2}; end; defval('N',128);
  % Steps of sizes that are being tried...
  if nargin>5; rindj=varargin{3}; end; defval('rindj',[2:2:N]);
  % Give this many experiments to each of the processors
  if nargin>6; Nspm=varargin{4}; end; defval('Nspm',5);
  % Axis handles to what will become the plot
  if nargin>7; ah=varargin{5}; end; defval('ah',gca);
  % Axis equalization parameters
  if nargin>8; gane=varargin{6}; end; 
  if nargin>9; miy=varargin{7}; end; 
  
  % Last default which wasn't an input
  defval('keepdata',0)
  
  % Initialize or accept the parameters and prefill some arrays
  [~,th0,p]=simulosl(th0,p);

  % For this set of parameters, make a unique hashed filename
  fname=hash([struct2array(orderfields(p)) th0 rindj Nspm abs(ptype)],'SHA-1');
  fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s.mat',fname));

  if ~exist(fnams,'file') 
    % Values and statistics that will be collected and kept
    [h,s,n,r,b,hm,hv,sm,nm,rm,h05,h95,s05,s95,n05,n95,r05,r95]=...
	deal(nan(length(rindj),1));
    thhat=nan(1,length(th0));
    NyNx=nan(length(rindj),2);
    if keepdata==1
      % Slots for one data example for everything
      H=cellnan([length(rindj) 1],rindj,rindj);
    end
    
    for index=1:length(rindj)
      % Change the size of the (square) patch under consideration
      NyNx(index,1:2)=rindj(index);
      p.NyNx=NyNx(index,:);

      tic
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(ptype,'poor')
        % Initialize to save time (?) I think I need COMPOSITE and not
        % CODISTRIBUTED since I do not need access to data between labs
        Hxv=Composite;
        for i=1:matlabpool('size'); Hxv{i}=nan(1,Nspm); end
	spmd
	  for sndex=1:Nspm
	    % Simulate new data with the same parameters and record the empirical variance
	    Hx=simulosl(th0,p);
	    % Keep track of the variance --- per processor, so these will be composites
	    Hxv(sndex)=var(Hx);
           end
        end
	% Flatten all the values since they sit uncomfortably in Composite class.
	try
	  Hxva=[Hxv{:}]; clear Hxv; 
	catch
	  Hxva=Hxv(:); clear Hxv; 
          end
	  % The stats of the poor-variances over the processors
	  hm(index)=mean(Hxva); hv(index)=var(Hxva);
	  h05(index)=prctile(Hxva,05); h95(index)=prctile(Hxva,95);
	  % Keep ONE poor-variance for every data patch size
	  h(index)=Hxva(1);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(ptype,'mle')
	% Initialize to save time (?) I think I need COMPOSITE and not
	% CODISTRIBUTED since I do not need access to data between labs
	thht=Composite;
	for i=1:matlabpool('size'); thht{i}=nan(Nspm,length(th0)); end
	% Now do a better job and run a real MLE inversion also
	spmd
	  for sndex=1:Nspm
	    % Simulate new data with the same parameters and record the
	    % EMPIRICAL variance of doing multiple simulations with these
	    % variables, rather than taking the word of "covh" for it
	    % even though we might want to take a peak at the output?
	    Hx=simulosl(th0,p);
	    try
	      % Make reasonbale guesses from the data themselves
	      [th,covh,~,~,scl]=mleosl(Hx,[var(Hx) 2.0 sqrt(prod(p.dydx.*p.NyNx))/5],p);
	      % It it was a NaN, fix the dimensions
	      if isnan(th); [th,scl]=deal(nan(1,length(th0))); end
	      % Output was scaled, so apply the scaling
	      thht(sndex,:)=th.*scl;
	    catch
	      thht(sndex,:)=nan(1,length(th0));
             end
           end
           end
	   % Flatten all the values since they sit uncomfortably in Composite class.
	   try
	     thhat=cat(1,thht{:}); clear thht; 
	   catch
	     thhat=thht; clear thht; 
           end      
	   % For the very small data sizes, how about some MLE cleanup?
	   % Because many times it's just not converging using FMINCON.
	   if rindj<15
	     thhat=trimit(thhat,80);
	     % Here we should save the number of survivors!
           end 
	   % The stats of the MLE-variances over the processors
	   sm(index)=nanmean(thhat(:,1));  sv(index)=nanvar(thhat(:,1));
	   nm(index)=nanmean(thhat(:,2));  nv(index)=nanvar(thhat(:,2));
	   rm(index)=nanmean(thhat(:,3));  rv(index)=nanvar(thhat(:,3));
	   s05(index)=prctile(thhat(:,1),05); s95(index)=prctile(thhat(:,1),95);
	   n05(index)=prctile(thhat(:,2),05); n95(index)=prctile(thhat(:,2),95);
	   r05(index)=prctile(thhat(:,3),05); r95(index)=prctile(thhat(:,3),95);
	   % Keep ONE MLE-variance for every data patch size
	   s(index)=thhat(1,1);
	   n(index)=thhat(1,2);
	   r(index)=thhat(1,3);
      end
      toc
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if keepdata==1
	% Keep one of the data sets in case you want to show it later
	H{index}=Hx{1}; clear Hx
      else
	H=NaN;
      end
      
      % Now predict the bias of the 'poor' variance from the known correlation structure
      b(index)= varbias(th0,p,1);
      b3(index)=varbias(th0,p,3);
    end

    % If we're worried, save inside the iteration?

    % Don't misinterpret the fact that we are saving a lot of different NyNx values
    rmfield(p,'NyNx');
    save(fnams,...
	 'h','s','n','r','b','b3','hm','hv','sm','nm','rm',...
	 'h05','h95','s05','s95','n05','n95','r05','r95',...
	 'p','th0','NyNx','thhat','Nspm','H')
  else
    disp(sprintf('Loading %s',fnams))
    load(fnams)
  end

  % Make the plots, in physical space (km!)
  struct2var(p)
  xlox=rindj*dydx(1)/1e3;
  xloy=linspace(0,[rindj(end)+1]*dydx(1),100)/1e3;
  xloy=unique([xloy linspace(0,dydx(1),100)/1e3]);
  
  % y-limits are based on the poor estimate
  defval('gane',range([h05 ; h ; h95])/20)
  defval('miy',[min([h ; h05]) max([h ; h95])]+[-gane gane]);
  
  hold on
  % Vertical and horizontal guides
  plot([th0(3) th0(3)]/1e3*pi*2,miy,'k-')
  %plot(xlim,[th0(1) th0(1)],'k--')

  switch ptype
   case 'poor'
    % Percentiles of the variances over all the realizations
    for index=1:length(rindj)
      pb(index)=plot([xlox(index) xlox(index)],[h05(index) h95(index)]);
    end
    % Rather plot a scaled version of the spatial-domain covariance itself!
    pp(4)=plot(xloy,maternosy(xloy*1e3,th0));
    % Variance of ONE of the realizations
    pp(1)=plot(xlox,h,'k');
    % Predicted mean of the variances knowing theoretical bias 
    pp(3)=plot([1*dydx(1)/1e3 xlox],th0(1)-[th0(1) ; b]);
    % Almost always awesome approximate prediction
    ppx=plot([1*dydx(1)/1e3 xlox],th0(1)-[th0(1) ; b3(:)],'kx-');
    % Mean of the variances over all the realizations
    pp(2)=plot(xlox,hm,'ko');
   case 'mle'
    % Percentiles of the variances over all the realizations
    for index=1:length(rindj)
      pb(index)=plot([xlox(index) xlox(index)],[s05(index) s95(index)]);
    end
    % Rather plot a scaled version of the spatial-domain covariance itself!
    pp(4)=plot(xloy,maternosy(xloy*1e3,th0));
    % MLE-variance of ONE of the realizations
    pp(1)=plot(xlox,s,'k');
    % Predicted mean of the variances knowing theoretical bias to be zero
    pp(3)=plot([1*dydx(1)/1e3 xlox],repmat(th0(1),1,length(xlox)+1));
    % Mean of the variances over all the realizations
    pp(2)=plot(xlox,sm,'ko');
  end

  hold off

  % Cosmetics
  set(pp(2),'Color','k','Marker','o','MarkerFaceC','w','MarkerEdgeC','k','MarkerS',4)
  set(pb,'Color',grey,'LineW',0.75)
  % Match the color in EGGERS1 with EGGERS4
  set(pp(3),'Color','b','LineW',1)
  set(pp(4),'Color','m','LineW',1)
  xlim([0 [rindj(end)+1]*dydx(1)]/1e3); 
  ylim(miy)
  longticks; 
  if length(ah)==1; shrink(gca); end
  t=title(sprintf('SIMULOSL with blurs = %i, var(Hx) versus %s^2',...
		  p.blurs,'\sigma'));
  movev(t,gane); grid off; 
  % xlabel('data size | distance (km)');
  xlabel(sprintf('grid size (km) ; 1 pixel = %i km',sqrt(prod(dydx))/1e3));
  %xl=xlabel(sprintf('grid size (km) ; %s = %i pixels',...
  %                  '\pi\rho',round(pi*th0(3)/sqrt(prod(dydx)))));
  % yl=ylabel('$s^2 \quad | \quad \mathcal{C}(r)$'); set(yl,'Interp','LaTex')
  yl=ylabel(sprintf('observed to predicted <s^2>/%s^2','\sigma'));
  % Scale the axis, it really doesn't matter
  % set(gca,'ytick',[0 th0(1) max(th0(1)+th0(1)/10,indeks(ylim,2))],'ytickl',{0, 1,' '})
  tix=4; tox=1.5;
  set(gca,'ytick',linspace(0,tox*th0(1),tix),'ytickl',{linspace(0,tox,tix)})

  % % Last minute fixin's in case we hadn't properly ordered to begin with
  delete(t)
  % for index=1:length(pb)
  %   bottom(pb(index),gca)
  % end
  % top(pp(2),gca)
  % bottom(pp(4),gca)
  
  % ax=xtraxis(gca,[],[],[],th0(1),'s2',[]); set(ax,'FontN','symbol')
  ax=xtraxis(gca,[],[],[],th0(1),[],[]); xl=xlim(ax);
  axt=text(xl(2)+range(xl)/30,th0(1),texlabel('sigma^2'),'FontS',get(gca,'FontS')+1);
  axis on 
  longticks(ax)
  
  % Delete the awesome prediction... assuming we verified it still works,
  % which we did in EGGERS1
  try ; delete(ppx) ; end
  
  % Figure out what all you will want to slap legends on, as in EGGERS1
  legs=[pp(2) pp(3) pb(1) pp(1) pp(4)];
  
  fig2print(gcf,'portrait')
  figna=figdisp([],sprintf('demo4_%s',ptype),[],1);
  system(sprintf('epstopdf %s.eps',figna));
  system(sprintf('rm -f %s.eps',figna));
  
  % Return the output if requested
  varns={gane,miy,ax,legs};
  varargout=varns(1:nargout);
elseif strcmp(th0,'demo5')
  % Shift the inputs so you can default them but also supply them
  % No further inputs needed, but if you have them, you keep them
  defval('params',[])
  defval('xver',[])
  % No inputs needed, but if you had them, you should use them in this order
  th0=params; p=xver;

  % Simulate
  [Hx,th0,p]=simulosl(th0,p);
  
  % How about some dorky and really terrible spatial-covariance estimation?
  ex=[0:p.NyNx(2)-1]*p.dydx(2);
  wy=[0:p.NyNx(1)-1]*p.dydx(1);
  dxxp=xxpdist(ex(:),wy(:));
  % [EX,WY]=meshgrid(ex,wy);
  % dxxp2=xxpdist([EX(:) WY(:)]);
  bins=linspace(min(dxxp(:)),max(dxxp(:))+1,20);
  bons=linspace(0,max(dxxp(:))+1,100);
  Ky=nan(1,length(bins)-1);
  for index=1:length(bins)-1
    % Identify pall the distances
    [i,j]=find(dxxp>=bins(index) & dxxp<bins(index+1));
    % Superpoor covariance estimate
    Ky(index)=mean([Hx(i)-mean(Hx(i))].*[Hx(j)-mean(Hx(j))]);
    % Variance is small but covariance is neglected and bias is
    % huge. Forget it... it's a simple illustration.
    dKy(index)=std([Hx(i)-mean(Hx(i))].*[Hx(j)-mean(Hx(j))])/length(i);
  end
  % Plot the results
  boks=[bins(1:end-1)+[bins(2)-bins(1)]/2]/1e3;
  plot(boks,Ky,'o'); hold on
  xlabel('distance [km]')
  ylabel('poor covariance estimate')
  % Blow up of these 'unreal' error bars
  fax=2;
  for ind=1:length(Ky)
    plot([boks(ind) boks(ind)],[Ky(ind)-fax*dKy(ind) Ky(ind)+fax*dKy(ind)],'-');
  end
  % Overley the theoretical value
  plot(bons/1e3,maternosy(bons,th0),'r-');
  % Plot in where you think the 2/3rd point of the covariance should be
  plot([th0(3) th0(3)]/1e3*pi*2,ylim,'k-')
  % Marke the zero line
  plot(xlim,[0 0],'k--')
  hold off
elseif strcmp(th0,'demo6')
  % Shift the inputs so you can default them but also supply them
  % No further inputs needed, but if you have them, you keep them
  defval('params',[])
  defval('xver',[])
  % No inputs needed, but if you had them, you should use them in this order
  th0=params; p=xver;

  % Just note that you need a large size for this to work well. Also,
  % there are boundary symmetry effects unless we "quarter"

  % Different, very labor-intensive tack to get to spatial covariance
  % Simulate N numbers of time
  N=102;
  [Hx,th0,p]=simulosl(th0,p);
  Hxx=deal(nan(N,size(v2s(Hx),2)));
  Hxxp=deal(nan(N,1));

  for index=1:N
    [Hx,th0]=simulosl(th0,p);
    % Demean...
    Hx=Hx-mean(Hx);
    % one random row
    wy=randi(size(v2s(Hx),1));
    % Collect the "data" at increasing separation
    Hx=v2s(Hx);
    % Pull out this random row
    Hxx(index,:)=Hx(wy,:);
    % Record the first point of this random row separately
    Hxxp(index)=Hx(wy,1);
    % Not sure why the means of these specfic positions or even fields
    % aren't zero themselves, that is probably due to the blurring which
    % leaves a non-zero wavenumber. But let's say it's the distances that
    % count. Hence the importance of demeaning the fields prior to
    % analysis, I guess. 
  end
  
  % Treat these random rows as 
  % Construct the spatial covariance, see if you get close
  Ky=nan(1,size(v2s(Hx),2));
  for index=1:size(v2s(Hx),2)
    Ky(index)=mean([Hxx(:,index)-mean(Hxx(:,index))].*[Hxxp-mean(Hxxp)]);
  end
  % Plot the results
  bons=[0:size(v2s(Hx),2)-1]*p.dydx(2);
  plot(bons/1e3,Ky,'-o'); hold on
  plot(bons/1e3,maternosy(bons,th0),'r-'); 
  plot([th0(3) th0(3)]/1e3*pi,ylim,'k--')
  xlim([0 size(v2s(Hx),2).*p.dydx(2)]/1e3)
  plot(xlim,[0 0],'k--')
  hold off
  xlabel('distance [km]')
  ylabel('poor covariance estimate')
  title(sprintf('blurs %i quart %i sims %i',p.blurs,p.quart,N))
  figna=figdisp([],sprintf('%i_%i_%i',p.blurs,p.quart,N),[],1);
  system(sprintf('epstopdf %s.eps',figna));
end

% Plotting routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tls,cbs,xcb,xa]=plotit(aha,dats,nm,stronk,strink,unid)
axes(aha)
imagesc(reshape(dats,nm)); 
axis image
limc=halverange(dats,95,NaN);
% Later, be more sophisticated than this
if limc(1)>-1
  limc=round(10*limc)/10;
else
  limc=round(limc);
end
set(aha,'clim',limc)
tls=title(stronk);
xa=xlabel(sprintf('mean %+6.3f ; stdev %6.3f %s',...
		  mean(dats(:)),std(dats(:)),unid));
cbs=colorbar('ver');
axes(cbs)
xcb=ylabel(sprintf(strink,unid));
set(cbs,'ylim',limc)
