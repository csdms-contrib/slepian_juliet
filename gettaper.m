function Tx=gettaper(params,tnam,attr)
% Tx=GETTAPER(PARAMS,TNAM,ATTR)
%
% Returns a taper for the grid specified by params for a rectangular taper or 
% a cosine taper on a rectangular grid or approximately on an irregular boundary
% (internal or external).
%
% INPUT:
%
% params    The grid parameters, including NyNx and dydx. If a
%           params.mask is provided and attr is two elements, the taper is 
%           applied to the mask
% tnam      Provide a string that refers to the taper style:
%           'rect'   - boxcar, Tx padded by attr*100%
%           'cosine' - a attr*100% cosine (Tukey) taper calculated as a 2D 
%                      extension of Walden & Percival 1993 eq 209
%           'slepian'- TODO: consider building a wrapper of SLEPIAN_FOXTROT
% attr      Attributes specific to the taper style; if params.mask is a field
%           and attr has two elements (scl of the main mask and percentage to
%           scale scl for the second mask), we will taper the requested mask
%           with attr(2)>0 tapering in and attr(2)<0 tapering out
%
% OUTPUT:
%
% Tx        Values of the taper on the full parameter grid
%
% SEE ALSO:
%
% maskit
%
% EXAMPLE:
%
% params=[];params.NyNx=[64 64];params.dydx=[1 1];
% Tx=getTaper(params,'cosine',0.2);
% getTaper('demo1')
% 
% % Taper a mask (irregular boundary)
% getTaper('demo2',[],[0.8 0.2])
% % Taper an anti-mask
% getTaper('demo2',[],[0.6 -0.4])
%
% % Study how well we can recover the Matern parameters using a unit taper if
% % a rectangular or cosine taper has been applied to a known field
% getTaper('demo3',[],0)
% % Include the taper in the analysis
% getTaper('demo3',[],1)
%
% Last modified by olwalbert-at-princeton.edu, 05/12/2025

if ~isstr(params)
  % The default will be an array of ones that is the size of the grid
  defval(tnam,'rect')
  % 2-dimensions
  d=2;
  NyNx=params.NyNx;
  dydx=params.dydx;
  if strcmp(tnam,'rect')
    % Percentage to pad grid or existing mask by (scalar between 0 and 1)
    defval('attr',0);
    if isfield(params,'mask')
      % We are requesting to chop a mask by attr*100%
      error('Option is not built yet; clear the mask field and try again.')
    else
      % We want to be able to calculate a mask from a taper the size of the
      % grid; by default, an MXN array of ones for attr==0
      Tx=ones(NyNx);
      % Apply zeros around the attr*100% of the boundary rather than padding as
      % we do not want to alter the size of the grid
      % This seems a bit tricky if we want to compensate for differences in
      % spacing; maybe this
      % adydx=(attr./dydx)./norm((attr./dydx),'fro');
      adydx=attr.*dydx;
      numrc=ceil(adydx.*NyNx/d^2);
      Tx(  [1:numrc(1) NyNx(1)-numrc(1)+1:NyNx(1)],:)=deal(0);
      Tx(:,[1:numrc(2) NyNx(2)-numrc(2)+1:NyNx(2)])  =deal(0);
      % Normalize according to Walden & Percival 1993 eq 208a
      % normfact=sqrt(sum(Tx.^2,'all'));
      % We are in 2-D, so probably more like
      % normfact=sum(Tx.^2,'all');
      % Let's avoid this here because we take care of it in the analysis
      normfact=1;
      Tx=Tx./normfact;
    end
  elseif strcmp(tnam,'cosine')
    defval('attr',0.2)
    if isfield(params,'mask') & numel(attr)==2
      % We are requesting to buffer a mask scaled to attr(1) with another mask
      % scaled by scl.*(1-attr(2)), so that attr(2) acts in a similar way to
      % attr in other cases; i.e., attr(2) is the proportion of the data to taper
      scl=attr(1);
      sclb=attr(2);
      % Interior mask
        % Find the mask vertices for the provided scaling, and for the scaling
        % provided in ATTR; attr should be less than 1 if we are buffering
        % outwards and greater than 1 if buffering inwards
        [~,I,~,~,~,cr]=maskit(randn(params.NyNx),params,scl);
        [~,Ib,~,~,~,crb]=maskit(randn(params.NyNx),params,scl.*(1-sclb));
        % Create denser arrays of polygon vertices
        xv=cr(:,1); yv=cr(:,2);
        xvb=crb(:,1); yvb=crb(:,2);
        % spline is not good enough for this: 
        % yv=spline(cr(1:end-1,1),cr(1:end-1,2),xv);
        dnsfct=5; % number of times to densify the polygon by averaging
        for dnd=1:dnsfct
          xv0=[xv(1:end-1) (xv(1:end-1)+xv(2:end))/2]';
          xv =[xv0(:); xv(end)];
          yv0=[yv(1:end-1) (yv(1:end-1)+yv(2:end))/2]';
          yv =[yv0(:); yv(end)];
          xvb0=[xvb(1:end-1) (xvb(1:end-1)+xvb(2:end))/2]';
          xvb =[xvb0(:); xvb(end)];
          yvb0=[yvb(1:end-1) (yvb(1:end-1)+yvb(2:end))/2]';
          yvb =[yvb0(:); yvb(end)];
        end
        cr=[xv yv];
        crb=[xvb yvb];
        % Find the mask elements for each segment of the mask in between cr and
        % crb, creating pieces
        [X,Y]=meshgrid(1:params.NyNx(2),1:params.NyNx(1));
        ctr=ceil(params.NyNx./2);
        bufm=zeros(params.NyNx);
        pin=zeros(params.NyNx);
        pon=zeros(params.NyNx);
        % Averaging over polygons defined by vertices nearby, not only nearest
        nr=3*dnsfct;
        for ind=1:size(cr,1)-nr
          vx=reshape([cr(ind:ind+nr,1) crb(ind:ind+nr,1)]',[],1);
          vy=reshape([cr(ind:ind+nr,2) crb(ind:ind+nr,2)]',[],1);
          % We must be very careful in defining the polygon's boundaries or
          % elements will be skipped
          [cvx,cvy]=poly2cw(vx,vy);
          idx = boundary(cvx,cvy);
          [inp,onp]=inpolygon(X,Y,cvx(idx),cvy(idx));
          inp=double(inp);
          % Keep track of which elements have been included
          pin=pin+inp;
          pon=pon+onp;
          % Assign weights from 0 to 1 based on distance from center
          numin=sum(inp(:));
          if numin>0
            [inpi,inpj]=find(inp>0);
            inpd=sqrt((inpi-ctr(1)).^2+(inpj-ctr(2)).^2);
            if sclb>0
              % 1s toward interior
              ninpd=(inpd-max(inpd(:)))./(min(inpd(:))-max(inpd(:)));
            else
              % 1s toward exterior
              ninpd=(inpd-min(inpd(:)))./(max(inpd(:))-min(inpd(:)));
            end
            ininpd=sub2ind(params.NyNx,inpi,inpj);
            bufm(ininpd)=bufm(ininpd)+1/2*(1-cos(pi*ninpd));
          end
        end
        bufm=bufm./pin;
        bufm(isnan(bufm))=deal(0);
        if sclb>0
          Tx=double(bufm.*~Ib.*I+Ib);
        else
          Tx=double(bufm.*~I.*Ib+~Ib);
        end
        % Sometimes elements are still missed or there are sharp edges; let's 
        % apply a little smoothing
        N=3;
        Tx2=conv2(Tx,ones(N)./N^2,'same');
        wt=1;
        Tx=(Tx+wt.*Tx2)./(1+wt);
        %Tx=conv2(Tx,[0 1 0; 1 1 1; 0 1 0]./4);
        if sclb>0
          Tx=I.*Tx;
        else
          Tx=~I.*Tx;
          % remove edge effects from smoothing
          Tx([1 end],:)=deal(1);
          Tx(:,[1 end])=deal(1);
        end
      % % If the mask were approximately round, this becomes an easier problem:
      % bbb=zeros(params.NyNx);
      % [bii,bjj]=find(abs(I-Ibufin));
      % bbd=sqrt((bii-ctr(1)).^2+(bjj-ctr(2)).^2);
      % nbd=(bbd-max(bbd(:)))./(min(bbd(:))-max(bbd(:)));
      % ibbs=sub2ind(params.NyNx,bii,bjj);
      % bbb(ibbs)=1/2*(1-cos(pi*nbd));
      % clf;imagesc(v2s(bbb,params))
    else
      % Extend Walden & Percival 1993 eq 209 for a attr x 100% cosine taper
      % to 2-D 
      % If we ignore dydx,
      ys=1:NyNx(1);
      xs=1:NyNx(2);
      ayax=floor(attr.*NyNx);
      % Define pieces in x and y
      pys=[{ys(ys<=ayax(1)/2)}...
           {ys(ys> ayax(1)/2 & ys<NyNx(1)+1-ayax(1)/2)}...
           {ys(ys>=NyNx(1)+1-ayax(1)/2)}];
      pxs=[{xs(xs<=ayax(2)/2)}...
           {xs(xs> ayax(2)/2 & xs<NyNx(2)+1-ayax(2)/2)}...
           {xs(xs>=NyNx(2)+1-ayax(2)/2)}];
      % Apply calculation piecewise
      tys{1}=1/2*(1-cos(2*pi*pys{1}./(ayax(1)+1)));
      txs{1}=1/2*(1-cos(2*pi*pxs{1}./(ayax(2)+1)));
      tys{2}=ones(size(pys{2})); 
      txs{2}=ones(size(pxs{2})); 
      tys{3}=1/2*(1-cos(2*pi*(NyNx(1)+1-pys{3})./(ayax(1)+1)));
      txs{3}=1/2*(1-cos(2*pi*(NyNx(2)+1-pxs{3})./(ayax(2)+1)));

      % Or, if we were to include dydx, we would want something like this,
      % ys=1:dydx(1):NyNx(1).*dydx(1);
      % xs=1:dydx(2):NyNx(2).*dydx(2);
      % ayax=floor(attr.*dydx.*NyNx);
      % ayax2=floor(attr.*dydx.*NyNx./2);
      % ayax=floor(attr.*NyNx);
      % ayax2=floor(attr.*NyNx./2);
      % % Split xs and ys into 3 pieces as cells
      % pys=[{ys(ys<=ayax2(1))}...
      %      {ys(ys> ayax2(1) & ys<NyNx(1)*dydx(1)+1-ayax2(1))}...
      %      {ys(ys>=NyNx(1)*dydx(1)+1-ayax2(1))}];
      % pxs=[{xs(xs<=ayax2(2))}...
      %      {xs(xs> ayax2(2) & xs<NyNx(2)*dydx(2)+1-ayax2(2))}...
      %      {xs(xs>=NyNx(2)*dydx(2)+1-ayax2(2))}];
      % % Apply calculation piecewise
      % tys{1}=1/2*(1-cos(2*pi*pys{1}./(ayax(1)+1)));
      % txs{1}=1/2*(1-cos(2*pi*pxs{1}./(ayax(2)+1)));
      % tys{2}=ones(size(pys{2})); 
      % txs{2}=ones(size(pxs{2})); 
      % tys{3}=1/2*(1-cos(2*pi*(NyNx(1)*dydx(1)+1-pys{3})./(ayax(1)+1)));
      % txs{3}=1/2*(1-cos(2*pi*(NyNx(2)*dydx(2)+1-pxs{3})./(ayax(2)+1)));

      % Collapse cell arrays
      tys=[tys{:}];
      txs=[txs{:}];
      % Expand to two-dimensions
      Tx=txs.*tys.';
      % Normalize according to Walden & Percival 1993 eq 208a
      normfact=sqrt(sum(Tx.^2,'all'));
      % BUT, we are in 2-dimensions, does the sqrt go away?
      normfact=sum(Tx.^2,'all');
      % BUTx2, we do this in the estimation procedure already. Ignore?
      normfact=1;
      % BUTx3 -- this still results in hat{sigma^2}/10 and hat{rho}/2; hat{nu}
      % is fine. Do we want to multiply by something instead so that it is as if
      % the grid sums to params.NyNx?
      normfact=sum(Tx,'all')/prod(params.NyNx);
      Tx=Tx./normfact;
    end
  elseif strcmp(tnam,'slepian')
    % Write a wrapper (eventually)
    error('We do not apply a slepian taper here; see SLEPIAN_FOXTROT instead')
  end
elseif strcmp(params,'demo1')
  % Visually inspect the 2D tapers and their profiles; compare to figures in 
  % Walden & Percival 1993 section 6.4, but remember that we have disabled
  % normalization
  p=[];p.NyNx=[64 64];p.dydx=[1 1];
  tnam='cosine';
  attr=0.2;
  Tx=getTaper(p,tnam,attr);

  clf
  [ah,ha,H]=krijetem(subnum(1,3));
  axes(ah(1))
  imagesc(Tx)
  axis image
  hold on
  yline(floor(size(Tx,1)/2),'w','LineWidth',2)
  xline(floor(size(Tx,2)/2),'w','LineWidth',2)
  xlabel('x2')
  ylabel('x1')
  longticks
  ti(1)=title(sprintf('2D %i%% %s taper',attr*100,tnam),'FontWeight','normal');

  axes(ah(2))
  plot(Tx(:,floor(end/2)),'ko')
  xline(0)
  axis square; axis padded
  ti(3)=title('1D profile along x1');
  xlabel('x1')
  ylabel('Tx')
  longticks
  xlim([0 p.NyNx(2)+1])

  axes(ah(3))
  plot(Tx(floor(end/2),:),'ko')
  xline(0)
  axis square; axis padded
  ti(2)=title('1D profile along x2');
  xlabel('x2')
  ylabel('Tx')
  longticks
  xlim([0 p.NyNx(1)+1])

  for ind=2:3
    movev(ti(ind),0.1*(ti(ind).Position(2)))
  end
  movev(ti(1),-6)
  sgtitle(sprintf('NyNx=[%i %i], dydx=[%0.2f %0.2f]',p.NyNx,p.dydx))
  saveas(gcf,sprintf('getTaper_demo1_%s',date),'epsc')
elseif strcmp(params,'demo2')
  % France example; we want this to request the smooth taper to be applied to a
  % mask by giving its name and attr. Let's hand attr as a vector, where
  % scl=attr(1) and attr=attr(2); if attr>0, then this is an interior mask, else
  % it is exterior
  defval('attr',[0.8 0.2]);
  scl=attr(1); 
  params=[];params.NyNx=[124 93];params.dydx=[1 1];params.mask='france';
  Tx=getTaper(params,'cosine',attr);
  [~,I,~,~,~,cr]  =maskit(randn(params.NyNx),params,scl);
  [~,Ib,~,~,~,crb]=maskit(randn(params.NyNx),params,scl.*(1-attr(2)));
  clf;imagesc(Tx);hold on;plot(cr(:,1),cr(:,2));plot(crb(:,1),crb(:,2));
  colorbar
  sgtitle(sprintf('NyNx=[%i %i], dydx=[%0.2f %0.2f], scl=%0.1f, attr=%0.1f',...
    params.NyNx,params.dydx,attr))
  saveas(gcf,sprintf('getTaper_demo2_%s',date),'epsc')
elseif strcmp(params,'demo3')
  % Study estimation of fields that have a simple taper applied with the option
  % to (attr==0) not include the taper in the estimation (i.e., unit taper) or to
  % (attr==1) include the true taper in the estimation
  defval(attr,0)
  % Let the taper grow and average over N simulations
  try
    load(sprintf('getTaper_demo3vars_%i.mat',attr))
  catch
    th=[10 1.5 5]; np=length(th);
    params=[]; params.NyNx=[128 128]; params.dydx=[1 1]; params.blurs=Inf; 
    if attr==0; params.taper=1; end
    M=8; N=10;
    thhatr=zeros(N,np,M);
    thhatc=zeros(N,np,M);
    for ind=1:M
      Tr=getTaper(params,'rect',0.1*ind);
      prcr(ind)=1-sum(Tr(:))./prod(size(Tr));
      Tc=getTaper(params,'cosine',0.05*ind);
      if attr==1; params.taper=Tc; end
      prcc(ind)=sum(Tc<max(Tc(:)),'all')./prod(size(Tc));
      for jnd=1:N
        Hx=simulosl(th,params);
        thhat=NaN;
        while isnan(thhat);
          [thhat,~,~,scl]=mleosl(Hx.*Tr(:),[],params,[],[],[],[1 1 1],1);
        end
        thhatr(jnd,:,ind)=thhat.*scl;
        thhat=NaN;
        while isnan(thhat);
          [thhat,~,~,scl]=mleosl(Hx.*Tc(:),[],params,[],[],th,[1 1 1],1);
        end
        thhatc(jnd,:,ind)=thhat.*scl;
      end
    end
    save(sprintf('getTaper_demo3vars_%i.mat',attr))
  end
  keyboard
  xlm=minmax([prcr.*100 prcc.*100])+[-2 2];
  clf;
  subplot(231)
  errorbar(prcr*100,squeeze(mean(thhatr(:,1,:),1)),squeeze(std(thhatr(:,1,:))./sqrt(N)))
  hold on; yline(th(1))
  ylabel('$\hat{\theta}$; rectangle taper','Interpreter','latex')
  title('$\sigma^2$','Interpreter','latex')
  ylim([1 11]); xlim(xlm)

  subplot(232)
  errorbar(prcr*100,squeeze(mean(thhatr(:,2,:),1)),squeeze(std(thhatr(:,2,:))./sqrt(N)))
  hold on; yline(th(2))
  title('$\nu$','Interpreter','latex')
  ylim([1.4 1.6]); xlim(xlm)

  subplot(233)
  errorbar(prcr.*100,squeeze(mean(thhatr(:,3,:),1)),squeeze(std(thhatr(:,3,:))./sqrt(N)))
  hold on; yline(th(3))
  title('$\rho$','Interpreter','latex')
  ylim([2 6]); xlim(xlm)

  subplot(234)
  errorbar(prcc.*100,squeeze(mean(thhatc(:,1,:),1)),squeeze(std(thhatc(:,1,:))./sqrt(N)))
  hold on; yline(th(1))
  ylabel('$\hat{\theta}$; cosine taper','Interpreter','latex')
  xlabel('percent tapered')
  ylim([1 11]); xlim(xlm)

  subplot(235)
  errorbar(prcc.*100,squeeze(mean(thhatc(:,2,:),1)),squeeze(std(thhatc(:,2,:))./sqrt(N)))
  hold on; yline(th(2))
  xlabel('percent tapered')
  ylim([1.4 1.6]); xlim(xlm)

  subplot(236)
  errorbar(prcc.*100,squeeze(mean(thhatc(:,3,:),1)),squeeze(std(thhatc(:,3,:))./sqrt(N)))
  hold on; yline(th(3))
  xlabel('percent tapered') 
  ylim([2 6]); xlim(xlm)
  sgtitle(sprintf('Estimating realizations of $\\theta_0=[%0.2f, %0.2f, %0.2f]$ with rectangle or cosine tapers applied from a unit taper',th),...
                  'Interpreter','latex')
  
end
