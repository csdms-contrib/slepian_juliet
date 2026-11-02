function varargout=maskit(v,p,scl,w,opt)
% [v,I,w,vw,scl,cr]=MASKIT(v,p,scl,w,opt)
%
% Arbitrary masking of a gridded spatial field.
%
% INPUT:
%
% v       A vector that is the unwrapping of a matrix
% p       A single parameter structure with AT LEAST this:
%           NyNx  number of samples in the y and x directions
%           mask an index matrix with a mask, OR:
%                 a region name that will be scaled to fit
%           mymx  a vector of indices in y and x that describe where to center
%                 the mask relative to the grid's center (null:  mymx=[0 0])
%           mr    angle (degrees) to rotate the mask by clockwise (null: mr=mod(360))
% scl     A scale between 0 and 1 with the fractional area that the
%         bounding box of the domain occupies; due to strict
%         inequality at maximum always leaves one pixel rim; and
%         not being used if an actual mask is input instead of a name
% w       A second vector that is the unwrapping of a second matrix
% opt     (Only for demo2) A field structure for optimization flags:
%           algo   the optimization algorithm
%           ifinv  ordered inversion flags for [s2 nu rho], e.g. [1 0 1] 
%
% OUTPUT:
%
% v       The first input with the samples inside I retained, rest set to NaN
% I       The mask, unrolled into a vector, note, this is an "anti-mask"
% w       The second input but everything outside the mask set to NaN
% vw      The merged field, second field w sucked into the masked first field v
% scl     A scale [0 to 1] expressing the non-missing data fraction
% cr      The colum,row index coordinates of the masking curve
%
% EXAMPLE:
%
% [~,~,I]=maskit('demo1'); % To get a quick mask to look at
% maskit('demo1') % A geographical region masking a single field
% maskit('demo1','england') % 'amazon', 'orinoco', for geographical variability
% maskit('demo2') % A geographical region merging two fields, with estimation
% maskit('demo2','england') % 'amazon', 'orinoco', for geographical variability
% maskit('demo3',[],rand) % Illustrating the nomenclature
% maskit('demo4') % Illustrates shifting and rotating by setting p.mymx and p.mr
% maskit('demo5') % Parameter recovery and uncertainty for (un)merged fields
%
% Last modified by olwalbert-princeton.edu, 09/08/2026
% Last modified by fjsimons-at-alum.mit.edu, 06/11/2024

% The default is the demo, for once
defval('v','demo1')
defstruct('opt',{'algo','ifinv'},{'unc',[1 1 1]});

if ~isstr(v)
    % Scale indicates survivorship (in symmetry with MUCKIT)
    % Default scale so bounding-box-area-proportionality is scl
    defval('scl',0.8)

    if (isstr(p.mask)|isfield(p,'XY'))
      % Get the mask...
      if isfield(p,'XY')
        % ... as XY vectors passed through the field p.XY
        XY=p.XY;
        % Scale relative to original XY
        cr(:,1)=             scale(p.XY(:,1),minmax(p.XY(:,1))+[1 -1]*(1-sqrt(scl))/2*(max(p.XY(:,1))-1));
        cr(:,2)=             scale(p.XY(:,2),minmax(p.XY(:,2))+[1 -1]*(1-sqrt(scl))/2*(max(p.XY(:,2))-1));
        %cr(:,2)=max(p.XY(:,2))-scale(p.XY(:,2),minmax(p.XY(:,2))+[1 -1]*(1-sqrt(scl))/2*(max(p.XY(:,2)-1))+1);
        % Check it some of the time
        if scl==1; diferm(p.XY-cr); end
      elseif isstr(p.mask)
        % ... by name 
        try 
          % first checking whether there is an associated *.m file (e.g., for 
          % Frederik's COAST masks) 
          XY=eval(p.mask);
        catch
          % and if not, generically load vectors from their *.mat file
          XY=getMask(p.mask);
        end
          % Scale the mask, remember for the curve the Y goes up, for the image, down
          % Isolated to the case of loading the mask by name, as we want to
          % preserve the option to NOT scale XY vectors relative to bounding box
          cr(:,1)=          scale(XY(:,1),[1 p.NyNx(2)]+[1 -1]*(1-sqrt(scl))/2*(p.NyNx(2)-1));
          cr(:,2)=p.NyNx(1)-scale(XY(:,2),[1 p.NyNx(1)]+[1 -1]*(1-sqrt(scl))/2*(p.NyNx(1)-1))+1;
      end
      % Optionally rotate the mask by p.mr, a clockwise rotation angle in degrees
      if isfield(p,'mr')
        mcr=mean(cr);
        cr=cr-mcr;
        R=[cos(p.mr*pi/180) -sin(p.mr*pi/180);...
           sin(p.mr*pi/180)  cos(p.mr*pi/180)];
        cr=(R*cr.').';
        cr=cr+mcr;
      end
      % Optionally shift the center of the mask relative to the center of the
      % grid
      if isfield(p,'mymx')
        % Shift positive vertical offset up and positive horizontal right in the 
        % axis ij convention
        % If we accidentally asked for a shift that is larger than the grid,
        % reset this value to the grid length
        p.mymx(p.mymx>p.NyNx) =p.NyNx(p.mymx>p.NyNx);
        p.mymx(p.mymx<-p.NyNx)=p.NyNx(p.mymx<-p.NyNx);
        % Shift the curve to shift the mask
        cr(:,1)=cr(:,1)+p.mymx(2);
        cr(:,2)=cr(:,2)-p.mymx(1);
      end
          
      % The underlying image grid
      [X,Y]=meshgrid(1:p.NyNx(2),1:p.NyNx(1));
          
      % Determine the mask
      I=inpolygon(X,Y,cr(:,1),cr(:,2));
    else
        % You already have it, and we won't make it any more complicated (skip
        % shift and rotation)
        cr=NaN;
        I=p.mask;
    end

    % Apply the mask to v, w and merge into vw as desired
    if nargin<4
        [w,vw]=deal(NaN);
    elseif nargout>3
        % Apply the mask to the second field
        w(~I(:))=NaN;
        vw=NaN;
        if nargout>4
            % Construct a new field with the SECOND field inside the mask and the FIRST outside
            vw=nan(size(v));
            vw(~~I(:))=w(~~I(:)); % avoid complaints about indices needing to
                                   %  be logicals in the case mask uses zeros
            vw(~I(:)) =v(~I(:));
        end
    end
    % Apply the mask to the first field
    v(~I(:))=NaN;
    % Note that you must not feed fields with NaN into MLEOSL - the taper I will
    % take of the selectivity
elseif strcmp(v,'demo1')
    % Capture the second input
    defval('p','france')
    defp=p; clear p
    % Now proceed with a fresh copy
    p.mask=defp;
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc; 
    [Hx,th,p]=simulosl([],p,1);
    % Capture the third input, default is none, which defaults inside
    defval('scl',[])
    [Hm,I,~,~,scl,cr]=maskit(Hx,p,scl);
    clf
    ah(1)=subplot(121); plotit(Hx,p,cr,th)
    ah(2)=subplot(122); plotit(Hm,p,cr,[])
    % Create right type of output even though there is non
    v=Hm; [w,vw]=deal(NaN);
elseif strcmp(v,'demo2')
    % Capture the second input
    defval('p','france')
    defp=p; clear p
    % Capture the fifth input
    defstruct('opt',{'algo','ifinv'},{[],[]});
    % Change the below if you like
    opt.ifinv=[1 1 1];
    
    % Now proceed with a fresh copy
    p.mask=defp;

    % Simulate using circulant embedding, no taper
    p.dydx=[1 1]; p.quart=0; p.blurs=Inf; p.kiso=NaN; clc;
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);
    % Something larger without overdoing it, check weirdness
    % If you want to generate the same hash you need to keep the same dimensions
    %    p.NyNx=[196 243];
    p.NyNx=[193 247];
    %p.NyNx=[193 236];

    % Do all the tests or not
    xver=0;
    % Decide on the length of the pause between displays, if any, 
    if xver==1
        pz=2;
    else
        pz=0.1;
        pz=1.5;
    end

    % Define the parameters for the two fields
    th1=[1.15 1.15 23000];
    th2=[2.90 2.75 11300];

    % Needing to be very explicit in order for parallel computing to work within the loop
    p.dydx=[10000 10000];
    p.taper=0; p.nugget=0;

    % Option to study the effect of setting kiso
    usekiso=0;
    if usekiso
      % value from eggers3
      p.kiso=pi./max(p.dydx);
    end

    % Option to use MLEOSL's default initialization
    usedini=1;
    if usedini
      thini1=[];
      thini2=[];
      aguess1=th1;
      aguess2=th2;
    end

    % Option to apply a smooth taper (~10% cosine for irregular boundary)
    usesmooth=0;
    % For some reason, cannot run in parfor for nonFrance; error hits in getMask

    % Number of processors, must agree with your machine
    NumWorkers=23;
    % Number of identical experiments to run experiments
    N=250;

    % Save some output for later? Make new filename hash from all relevant input
    fname=hash([struct2array(orderfields(p)) struct2array(orderfields(opt)) th1 th2 N],'SHA-1');
    
    % You need to have an environmental variable file structure set up
    fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s_%s.mat','MASKIT',fname));

    % Let's shift the mask about by setting offset parameter mymx at random
    p.mymx=randi(floor(sqrt(sum(p.NyNx.^2))/8),[1 2]);
    % Randomly select direction of shift
    p.mymx=p.mymx.*sign(randn(1,2));
    % Spin it, also at random
    p.mr=randi(360)*sign(randn(1));

    p.mymx=0; p.mr=0;

    % How much should the masked area occupy?
    %scl=[];
    % Generate an initial mask and smoothing taper outside of the loop 
    [~,I,~,~,scl]=maskit(randn(p.NyNx),p);
    % Whole field
    p1=p; p1.taper=0;
    if usesmooth; p1.mask=[]; Tx=gettaper(p1,'cosine',0.1); p1.taper=Tx; end
    % Masked field
    p2=p; p2.taper=~I;
    if usesmooth; Tx=gettaper(p2,'cosine',[scl 0.1]); p2.taper=Tx; p2.mask=Tx; end
    % Anti-masked field
    p3=p; p3.taper=I;
    if usesmooth; Tx=gettaper(p3,'cosine',[scl -0.1]); p3.taper=Tx; p3.mask=Tx; end
    % Mixed field
    p4=p; p4.taper=0;
    if usesmooth; p4.mask=[]; Tx=gettaper(p4,'cosine',0.1); p4.taper=Tx; end

    % Run the experiment anew
    if ~exist(fnams,'file') % || 1==1
        % Initialize the pool of workers
        if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end

        % Should we preallocate?
        
        % Run the experiment!
        %parfor index=1:N
        for index=1:N
            clc; disp(sprintf('\n Simulating all fields \n'))

            % Simulate first field
            Hx=simulosl(th1,p,1);
            % Simulate second field
            Gx=simulosl(th2,p,1);

            % Now do the masking and the merging
            [Hm,I,Gm,HG,scl,cr]=maskit(Hx,p,[],Gx);
disp('266')
            % So HG is the mixed field
            v=Hm; w=Gm; vw=HG;

            % As appropriate you'll force the use of BLUROSY in MATERNOSP in LOGLIOSL in MLEOSL
            % This is confusing PARFOR but works with FOR
            % p.blurs=-1;

            if ~usedini
              % Make a close initial guess?
              thini1i(index,:)=th1+(-1).^randi(2,[1 3]).*th1/1000;
              thini2i(index,:)=th2+(-1).^randi(2,[1 3]).*th2/1000;
              thini1=thini1i(index,:);
              thini2=thini2i(index,:);
              aguess1=[];
              aguess2=[];
            end

            % Don't perturb the thini for the one you don't want to move; could
            % have also have stuck in th1 and th2 into the aguess slot of
            % MLEOSL, whence they would be appropriately (or not) perturbed
            if opt.ifinv==[1 0 1]
                % This didn't work
                %thini1(index,2)=th1(2);
                %thini2(index,2)=th2(2);
                % But this did
                thini1(index,:)=th1+(-1).^randi(2,[1 3]).*[th1(1) 0 th1(3)]/1000;
                thini2(index,:)=th2+(-1).^randi(2,[1 3]).*[th2(1) 0 th2(3)]/1000;
            end
            
            % Recover the parameters of the full original fields without any masking
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            pause(pz); clc; disp(sprintf('\n Estimating first whole field \n'))
            [thhat4(index,:),~,~,scl4(index,:)]=mleosl(Hx,thini1,p1,opt.algo,[],aguess1,opt.ifinv,xver);
            pause(pz); clc ; disp(sprintf('\n Estimating second whole field \n'))
            [thhat5(index,:),~,~,scl5(index,:)]=mleosl(Gx,thini2,p1,opt.algo,[],aguess2,opt.ifinv,xver);

            % Now recover the parameters of the mixed field but only in the region ~I
            % Perform the optimization on the complement which should look like the FIRST field
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            pause(pz); clc; disp(sprintf('\n Estimating first partial field \n'))
            [thhat1(index,:),~,~,scl1(index,:)]=mleosl(HG,thini1,p2,opt.algo,[],aguess1,opt.ifinv,xver);
            
            % Now recover the parameters of the insert which should look like the SECOND field
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            pause(pz); clc; disp(sprintf('\n Estimating second partial field \n'))
            [thhat2(index,:),~,~,scl2(index,:)]=mleosl(HG,thini2,p3,opt.algo,[],aguess2,opt.ifinv,xver);

            % Now recover the parameters of the merged field without knowing of the partition
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            pause(pz); clc; disp(sprintf('\n Estimating whole mixed field \n'))
            % Force the inversion of all parameters in this particular case 
            [thhat3(index,:),~,~,scl3(index,:)]=mleosl(HG,[],             p4,opt.algo,[],[],[1 1 1],xver);

            % Pause so you can watch live (big pause if xver==0)
            pause(pz)

            % This is unnecessary for PARFOR but works with FOR
            % Reset to the original values
            %p.blurs=Inf; p.taper=0;
        end
        % Save all the output that you'll need later, may add more later
        save(fnams,'th1','th2','p','thhat1','thhat2','scl1','scl2',...
             'thhat3','scl3','thhat4','scl4','thhat5','scl5','thini1','thini2','N')
    else
        % Load all the saved output
        load(fnams)
    end

    % And now look at the statistics of the recovery on screen
    str0='%8.2f %5.2f %6.0f';
    str1=sprintf('%s %s\n',str0,str0);
    str2=sprintf('%s\n',str0);
    disp(sprintf('\nFirst partial | Second partial\n'))
    disp(sprintf(str1,[thhat1.*scl1 thhat2.*scl2]'))
    disp(sprintf('\nFirst whole | Second whole\n'))
    disp(sprintf(str1,[thhat4.*scl4 thhat5.*scl5]'))
    disp(sprintf('\nMixed whole\n'))
    disp(sprintf(str2,[thhat3.*scl3]'))
    
    % Save the figures bjootifooly
    % Remake one sampled field just to make the visual
    [Hx,~,~,~,~,Sb1]=simulosl(th1,p,1);
    [Gx,~,~,~,~,Sb2]=simulosl(th2,p,1);
    [~,~,~,HG,scl,cr]=maskit(Hx,p,[],Gx);

    % This was confusing PARFOR so we took it out of the loop
    if xver==1
        % Very explicit comparison like BLUROSY
        p.blurs=-1; Sbar1=blurosy(th1,p,1); 
        p.blurs=-1; Sbar2=blurosy(th2,p,1);
        % Are any of them real bad?
        if max(Sb1./Sbar1)>20; warning('You will want to review this code') ; end
        if max(Sb2./Sbar2)>20; warning('You will want to review this code') ; end
        % Reset to the original value
        p.blurs=Inf;
    end

    % Make a visual for good measure
    figure(3)
    clf
    if usesmooth; Hx=Hx.*Tx(:); Gx=Gx.*Tx(:); end
    ah(1)=subplot(221); plotit(Hx,p,cr,th1)
    ah(1).Title.Position(2)=ah(1).Title.Position(2)+7;
    ah(3)=subplot(223); plotit(Gx,p,cr,th2)
    ah(3).Title.Position(2)=ah(3).Title.Position(2)+7;
    ah(4)=subplot(224); plotit(HG,p,cr,[])
    movev(ah(4),0.25)
    figdisp(sprintf('%s_1',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)

    % Then plot the estimates using MLEPLOS, remember second field is the inserted one
    % So this is the ANTI region
    figure(1)
    mleplos(thhat1.*scl1,th1,[],[],[],[],[],p2,sprintf('anti %s',p.mask),[],opt.ifinv)
    pause(pz)
    figure(1)
    figdisp(sprintf('%s_2a',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    pause(pz)
    clf
    figure(2)
    pause(pz)
    figdisp(sprintf('%s_2b',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    pause(pz)
    clf
    % And this is the SELECTED region
    figure(1)
    mleplos(thhat2.*scl2,th2,[],[],[],[],[],p3,p.mask,[],opt.ifinv)
    pause(pz)
    figure(1)
    pause(pz)
    figdisp(sprintf('%s_3a',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    pause(pz)
    clf
    figure(2)
    pause(pz)
    figdisp(sprintf('%s_3b',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    pause(pz)

    ptcldfig=0;
    if ptcldfig
      % examine estimates in a 3D scatter plot
      figure(4)
      clf
      tss={thhat1.*scl1;thhat2.*scl2;thhat3.*scl3;thhat4.*scl4;thhat5.*scl5};
      css={'#0072BD','#0072BD','#77AC30','#EDB120','#EDB120'};
      mss={'^','o','o','^','o'};
      hold on
      scatter3(th1(:,1),th1(:,2),th1(:,3),1000,'kx')
      scatter3(th2(:,1),th2(:,2),th2(:,3),1000,'k*')
      for tnd=1:5
        ts=tss{tnd};
        scatter3(ts(:,1),ts(:,2),ts(:,3),...
                 mss{tnd},'MarkerFaceColor',css{tnd},'MarkerEdgeColor',css{tnd},...
                 'MarkerFaceAlpha',0.6)                                   
      end
      view([-18 40])
      grid on
      xlabel('\sigma^2')                                                              
      ylabel('\nu')                                                                   
      zlabel('\rho')     
      legend({'th1','th2','first partial (anti)','second partial (mask)',...
              'mixed whole','first whole','second whole'})
      title(sprintf('%i MLE Simulations; %ix%i grid; %ix%i km; blur %s; %s',...
            size(thhat1,1),p.NyNx,p.NyNx.*p.dydx./1000,p.blurs,upper(p.mask)))
      figdisp(sprintf('%s_4',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    end

    % With PARFOR none of the once-used are available out of the loop but in
    % this demo2 you don't want any output anyway, so put in empties
    [v,I,w,vw,scl,cr]=deal(NaN);
elseif strcmp(v,'demo3')
    th1=[1.15 1.15 23000];
    th2=[2.90 2.75 11300];
    p.mask='france';
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);
    % Capture the third input, default is none, which defaults inside
    defval('scl',[])
    
    % Remake one sampled field just to make the visual
    [Hx,~,p,~,~,Sb1]=simulosl(th1,p);
    [Gx,~,~,~,~,Sb2]=simulosl(th2,p);
    [v,I,w,vw,scl,cr]=maskit(Hx,p,scl,Gx);

    % Make a visual for good measure
    clf
    ah(1)=subplot(221); plotit(v,p,cr,sprintf('%i %% v',round(100*scl)))
    ah(3)=subplot(223); plotit(w,p,cr,sprintf('%i %% w',round(100*scl)))
    ah(4)=subplot(224); plotit(vw,p,cr,'vw')
    movev(ah(4),0.25)
elseif strcmp(v,'demo4')
    % Demonstration of setting p.mymx and p.mr for MASKIT: 
    % positive p.mymx(1) translates the mask up, positive p.mymx(2) translates 
    % mask to the right in the axis  ij convention
    % p.mr is the clockwise angle to rotate, in degrees
    % Select grid parameters at random
    p=[];p.NyNx=randi(120,[1,2])+36;p.dydx=[1 1];p.blurs=-1;
    p.mask='france';
    mscl=(randi(80)+20)/100;

    % Create the mask without mymx
    [v,I,w,vw,scl,cr]=maskit(rand(p.NyNx),p,mscl);

    % Set mymx at random
    p.mymx=randi(40,[1 2]);
    % Randomly select direction of shift
    p.mymx=p.mymx.*sign(randn(1,2));
    % Create the mask with mymx
    [~,Imymx]=maskit(rand(p.NyNx),p,mscl);

    % Give it a spin
    p.mr=randi(360)*sign(randn(1));
    [~,Imymxmr]=maskit(rand(p.NyNx),p,mscl);

    % Plot the original grid-centered mask, the shifted mask, and the shifted,
    % rotated mask
    clf;
    subplot(131);imagesc(I)
    subplot(132);imagesc(Imymx)
    title(sprintf('vertical: %i; horizontal: %i',p.mymx))
    subplot(133);imagesc(Imymxmr)
    title(sprintf('rotation: %i%s',p.mr,char(176)))

    % Alternatively, if we wanted to allow the mask to wrap around,
    % mymx=mod(floor(p.mymx),p.NyNx./2);
    % I=circshift(I,p.mymx.*[-1 1]);
elseif strcmp(v,'demo5')
    % Demonstrate recovered parameters as histograms for a variety of mask
    % options: (a) full field I / II, (b) the mixed field first with I interior
    % then II interior, (c) tapered field with I / II interior, (d) tapered 
    % field with I / II exterior (the anti-mask), (e&f) same as (c) and (d) but
    % analyzed with a unit taper, (g&h) same as (b) and (c) but smooth taper 
    % applied to boundary of mask and then analyzed with a square taper
    % msk='france'; datum='08-May-2025';%'30-Apr-2025'; 
    msk='newjersey'; datum='02-Jun-2025';
    msk='newjersey'; datum='25-Aug-2025';%'24-Aug-2025';
    % msk='doves'; datum='12-May-2025';
    np=3;
    % We will need to run several experiments; let's save them for future
    % loading
    try
      % Load the experiments if we have alread run them
      load(sprintf('maskit_demo5vars_%s_%s.mat',msk,datum));
      keyboard
      % a few more variables that will be needed in the figure 
      numreals=size(thai,1);
      simparams=params;simparams.blurs=Inf;
      % the full fields
      Hai =simulosl(th1,simparams,1);
      Haii=simulosl(th2,simparams,1);
      % the mixed fields
      scl=[];
      [Hci,I,Hcii,Hbi,scl,cr]=maskit(Hai,simparams,scl,Haii);
      [~,~,~,Hbii,~,~]=maskit(Haii,simparams,scl,Hai);
      % the field tapered by the mask
      Hci=Hai.*I(:);
      Hcii=Haii.*I(:);
      % the field tapered by the anti-mask
      Hdi =Hai.*~I(:);
      Hdii=Haii.*~I(:);
      % the smoothly tapered mask
      Txm=gettaper(params,'cosine',[scl 0.2]); 
      % the smoothly tapered anti-mask
      Txa=gettaper(params,'cosine',[scl -0.2]);
      
      % the smoothly tapered fields, in space
      Hgi =Hai.*Txm(:);
      Hgii=Haii.*Txm(:);
      Hhi =Hai.*Txa(:);
      Hhii=Haii.*Txa(:);
    catch
      % Run the experiments in a similar fashion as DEMO2
      numreals=200;
      th1=[1.15 0.65 1800];
      th2=[2.90 1.75 3000];
      params=[];params.NyNx=[124 93];params.dydx=[1e3 1e3]; params.blurs=-1;
      params.mask=msk;
      scl=[];
      simparams=params;simparams.blurs=Inf;
      opt.algo='unc';opt.ifinv=[1 1 1];xver=0;
      % Calculate tapers that will be applied either directly to the data or in
      % MLEOSL
      % A unit-taper (square):    pu
      % A taper for the mask:     pm
      % A taper for the antimask: pa
      % Smooth taper to be directly applied to the data (mask): Txm
      % Smooth taper to be directly applied to the data (anti): Txa
      [~,I,~,~,scl]=maskit(randn(simparams.NyNx),simparams,scl);
      pu=params;pu.taper= 1;
      pm=params;pm.taper= double(I);
      pa=params;pa.taper=double(~I);
      % IM=double(I);  IM(IM==0)=deal(NaN); IM=IM(:);
      % IA=double(~I); IA(IA==0)=deal(NaN); IA=IA(:);
      Txm=gettaper(pm,'cosine',[scl 0.2]); 
      Txa=gettaper(pa,'cosine',[scl -0.2]);
      psm=params;psm.taper=Txm;
      psa=params;psa.taper=Txa;

      for index=1:numreals
        clc; disp(sprintf('\n Simulating all fields \n'))
        % Simulate full field I
        Hai =simulosl(th1,simparams,1);
        % Simulate full field II
        Haii=simulosl(th2,simparams,1);
              
        % Now do the masking and the merging
        % Full field I: Hai; full field II: Haii; 
        % Merged field I inside: Hbi; merged field II inside: Hbii
        % Masked field I: Hci; masked field II: Hcii
        % Anti-masked field I: Hdi; anti-masked field II: Hdii
        % Smoothly tapered mask I: Hgi; smoothly tapered mask II: Hgii
        % Smoothly tapered anti-mask I: Hhi; smoothly tapered anti-mask II: Hhii
        [Hci,I,Hcii,Hbi,scl,cr]=maskit(Hai,simparams,scl,Haii);
        [~,~,~,Hbii,~,~]=maskit(Haii,simparams,scl,Hai);
        Hci=Hai.*I(:);
        Hcii=Haii.*I(:);
        Hdi =Hai.*~I(:);
        Hdii=Haii.*~I(:);

        %%%
        % % Apply the smoothing taper to mask buffer
        % attr=0.01;Ns=sqrt(sum(params.NyNx.^2));
        % buf=attr.*Ns/2; % buf is measured in arc distance degrees
        % npts=25; % 13 points are used by default for each vertex; make this smoother
        % xsp=cr(:,1); ysp=cr(:,2);
        % % spline is not good: ysp=spline(cr(1:end-1,1),cr(1:end-1,2),xsp);
        % dnsfctn=4; % number of times to densify the polygon by averaging
        % for dnd=1:dnsfctn
        %   xsp0=[xsp(1:end-1) (xsp(1:end-1)+xsp(2:end))/2]';
        %   xsp =[xsp0(:); xsp(end)];
        %   ysp0=[ysp(1:end-1) (ysp(1:end-1)+ysp(2:end))/2]';
        %   ysp =[ysp0(:); ysp(end)];
        % end
        % [xb,yb]=bufferm(xsp,ysp,abs(buf),'in',npts);
        % % Bufferm selects an uneven amount of cells in different regions of the
        % % polygon. I tried making the mask denser by averaging between points 
        % % (splines was not so nice), but this wasn't sufficient to make the
        % % buffer zone even, nor was increasing the number of points at each
        % % vertex (npts)
        % mcr=[xb yb];
        % % The underlying image grid
        % [X,Y]=meshgrid(1:params.NyNx(2),1:params.NyNx(1));
        % % Determine the mask
        % tmsmooth=inpolygon(X,Y,mcr(:,1),mcr(:,2));
        % % imagesc(int16(inpolygon(X,Y,cr(:,1),cr(:,2)))+2.*int16(tmsmooth))
        % % For the antimask
        % [yb,xb]=bufferm(cr(:,2),cr(:,1),abs(buf),'in',npts);
        % acr=[xb yb];
        % % Determine the mask
        % inp=inpolygon(X,Y,acr(:,1),acr(:,2));
        % % The polygon does not step in evenly.
        %%%

        % Instead, let's apply the taper calculated from nested masks
        Hgi =Hai.*Txm(:);
        Hgii=Haii.*Txm(:);
        Hhi =Hai.*Txa(:);
        Hhii=Haii.*Txa(:);

        % Make a close initial guess?
        % thini1=th1+(-1).^randi(2,[1 3]).*th1/1000;
        % thini2=th2+(-1).^randi(2,[1 3]).*th2/1000;
        thini1=[];
        thini2=[];

        % Experiment A: Full fields
        % Recover the parameters of the full original fields without any masking
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hai ,thini1,pu,opt.algo,[],th1,opt.ifinv,xver);
        end
        thai(index,:) =thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Haii,thini2,pu,opt.algo,[],th2,opt.ifinv,xver);
        end
        thaii(index,:)=thhat.*mcl;

        % Experiment B: Merged fields
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hbi ,[],pu,opt.algo,[],th1,[1 1 1],xver);
        end
        thbi(index,:) =thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hbii,[],pu,opt.algo,[],th2,[1 1 1],xver);
        end
        thbii(index,:)=thhat.*mcl;

        % Experiment C: Masked fields
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          % We don't want the taper applied twice inside of mleosl, something
          % funny happens and we can modify this later
          %[thhat,~,~,mcl]=mleosl(Hci ,thini1,pm,opt.algo,[],th1,opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Hai ,thini1,pm,opt.algo,[],th1,opt.ifinv,xver);
        end
        thci (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          % We don't want the taper applied twice inside of mleosl, something
          % funny happens and we can modify this later; instead, provide full
          % field and selective taper that will be used
          %[thhat,~,~,mcl]=mleosl(Hcii,thini2,pm,opt.algo,[],th2,opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Haii,thini2,pm,opt.algo,[],th2,opt.ifinv,xver);
        end
        thcii(index,:)=thhat.*mcl;
          
        % Experiment D: Anti-masked fields
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          %[thhat,~,~,mcl]=mleosl(Hdi ,thini1,pa,opt.algo,[],[],opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Hai ,thini1,pa,opt.algo,[],[],opt.ifinv,xver);
        end
        thdi (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          %[thhat,~,~,mcl]=mleosl(Hdii,thini2,pa,opt.algo,[],[],opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Haii,thini2,pa,opt.algo,[],[],opt.ifinv,xver);
        end
        thdii(index,:)=thhat.*mcl;

        % Experiment E: Masked fields, but estimated on a grid with a unit taper 
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hci ,thini1,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thei (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hcii,thini2,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        theii(index,:)=thhat.*mcl;

        % Experiment F: Anti-masked fields, but estimated on a grid with a unit 
        % taper 
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hdi ,thini1,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thfi (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hdii,thini2,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thfii(index,:)=thhat.*mcl;

        % Experiment G: Masked fields, smoothly tapered, but estimated on a grid
        % with a unit taper 
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hgi ,thini1,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thgi (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hgii,thini2,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thgii(index,:)=thhat.*mcl;

        % Experiment H: Anti-masked fields, smoothly tapered, but estimated on a
        % grid with a unit taper 
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hhi ,thini1,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thhi (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          [thhat,~,~,mcl]=mleosl(Hhii,thini2,pu,opt.algo,[],[],opt.ifinv,xver);
        end
        thhii(index,:)=thhat.*mcl;
        
        % Experiment I: Masked field estimated with a smoothed mask taper
        % remember that the smooth taper is applied to the field when converting
        % to the spectral domain
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          %[thhat,~,~,mcl]=mleosl(Hci ,thini1,pm,opt.algo,[],[],opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Hai ,thini1,psm,opt.algo,[],[],opt.ifinv,xver);
        end
        thii (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          %[thhat,~,~,mcl]=mleosl(Hcii,thini2,pm,opt.algo,[],[],opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Haii,thini2,psm,opt.algo,[],[],opt.ifinv,xver);
        end
        thiii(index,:)=thhat.*mcl;

        % Experiment J: Anti-masked fields, smoothly tapered, estimated with the 
        % anti-mask taper
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          %[thhat,~,~,mcl]=mleosl(Hdi ,thini1,pa,opt.algo,[],[],opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Hai ,thini1,psa,opt.algo,[],[],opt.ifinv,xver);
        end
        thji (index,:)=thhat.*mcl;
        %
        thhat=NaN;mcl=NaN;
        while isnan(thhat.*mcl)
          %[thhat,~,~,mcl]=mleosl(Hdii,thini2,pa,opt.algo,[],[],opt.ifinv,xver);
          [thhat,~,~,mcl]=mleosl(Haii,thini2,psa,opt.algo,[],[],opt.ifinv,xver);
        end
        thjii(index,:)=thhat.*mcl;
      end

      % Calculate the empirical covariance from each experiment
      cvaie=nancov(thai); cvaiie=nancov(thaii); cvbie=nancov(thbi); cvbiie=nancov(thbii);
      cvcie=nancov(thci); cvciie=nancov(thcii); cvdie=nancov(thdi); cvdiie=nancov(thdii);
      cveie=nancov(thei); cveiie=nancov(theii); cvfie=nancov(thfi); cvfiie=nancov(thfii);
      cvgie=nancov(thgi); cvgiie=nancov(thgii); cvhie=nancov(thhi); cvhiie=nancov(thhii);
      cviie=nancov(thii); cviiie=nancov(thiii); cvjie=nancov(thji); cvjiie=nancov(thjii);
      % Calculate the analytical covariance for each Matern parameter and taper
      % pair  
      cm=2;
      cvaid =covthosl(mean(thai),pu,cm,opt.ifinv);  
      cvaiid =covthosl(mean(thaii),pu,cm,opt.ifinv);
      cvbid =covthosl(mean(thbi),pu,cm,opt.ifinv);  
      cvbiid =covthosl(mean(thbii),pu,cm,opt.ifinv);
      cvcid =covthosl(mean(thci),pm,cm,opt.ifinv);  
      cvciid =covthosl(mean(thcii),pm,cm,opt.ifinv);
      cvdid =covthosl(mean(thdi),pa,cm,opt.ifinv);  
      cvdiid =covthosl(mean(thdii),pa,cm,opt.ifinv);
      cveid =covthosl(mean(thei),pu,cm,opt.ifinv);  
      cveiid =covthosl(mean(theii),pu,cm,opt.ifinv);
      cvfid =covthosl(mean(thfi),pu,cm,opt.ifinv);  
      cvfiid =covthosl(mean(thfii),pu,cm,opt.ifinv);
      cvgid =covthosl(mean(thgi),pu,cm,opt.ifinv);  
      cvgiid =covthosl(mean(thgii),pu,cm,opt.ifinv);
      cvhid =covthosl(mean(thhi),pu,cm,opt.ifinv);  
      cvhiid =covthosl(mean(thhii),pu,cm,opt.ifinv);
      cviid =covthosl(mean(thii),psm,cm,opt.ifinv);  
      cviiid =covthosl(mean(thiii),psm,cm,opt.ifinv);
      cvjid =covthosl(mean(thji),psa,cm,opt.ifinv);  
      cvjiid =covthosl(mean(thjii),psa,cm,opt.ifinv);

      % Save the results
      save(sprintf('maskit_demo5vars_%s_%s',msk,datum),...
          'params','th1','th2',...
          'thai','thaii','cvaie','cvaiie','cvaid','cvaiid',...
          'thbi','thbii','cvbie','cvbiie','cvbid','cvbiid',...
          'thci','thcii','cvcie','cvciie','cvcid','cvciid',...
          'thdi','thdii','cvdie','cvdiie','cvdid','cvdiid',...
          'thei','theii','cveie','cveiie','cveid','cveiid',...
          'thfi','thfii','cvfie','cvfiie','cvfid','cvfiid',...
          'thgi','thgii','cvgie','cvgiie','cvgid','cvgiid',...
          'thhi','thhii','cvhie','cvhiie','cvhid','cvhiid',...
          'thii','thiii','cviie','cviiie','cviid','cviiid',...
          'thji','thjii','cvjie','cvjiie','cvjid','cvjiid');
    end
    outrm=1;
    if outrm
      thai=madout(thai,1); thaii=madout(thaii,1);
      thbi=madout(thbi,1); thbii=madout(thbii,1);
      thci=madout(thci,1); thcii=madout(thcii,1);
      thdi=madout(thdi,1); thdii=madout(thdii,1);
      thei=madout(thei,1); theii=madout(theii,1);
      thfi=madout(thfi,1); thfii=madout(thfii,1);
      thgi=madout(thgi,1); thgii=madout(thgii,1);
      thhi=madout(thhi,1); thhii=madout(thhii,1);
      thii=madout(thii,1); thiii=madout(thiii,1);
      thji=madout(thji,1); thjii=madout(thjii,1);
      % recalculate covariances
      % Calculate the empirical covariance from each experiment                 
      cvaie=nancov(thai); cvaiie=nancov(thaii); cvbie=nancov(thbi); cvbiie=nancov(thbii);
      cvcie=nancov(thci); cvciie=nancov(thcii); cvdie=nancov(thdi); cvdiie=nancov(thdii);
      cveie=nancov(thei); cveiie=nancov(theii); cvfie=nancov(thfi); cvfiie=nancov(thfii);
      cvgie=nancov(thgi); cvgiie=nancov(thgii); cvhie=nancov(thhi); cvhiie=nancov(thhii);
      cviie=nancov(thii); cviiie=nancov(thiii); cvjie=nancov(thji); cvjiie=nancov(thjii);
      % Calculate the analytical covariance for each Matern parameter and taper 
      % pair                                                                    
      pu=params;pu.taper= 1;
      pm=params;pm.taper= double(I);
      pa=params;pa.taper=double(~I);
      Txm=gettaper(pm,'cosine',[scl 0.2]); 
      Txa=gettaper(pa,'cosine',[scl -0.2]);
      psm=params;psm.taper=Txm;
      psa=params;psa.taper=Txa;
      cm=2;                                                                     
      cvaid =covthosl(mean(thai),pu,cm,opt.ifinv);                              
      cvaiid =covthosl(mean(thaii),pu,cm,opt.ifinv);                            
      cvbid =covthosl(mean(thbi),pu,cm,opt.ifinv);                              
      cvbiid =covthosl(mean(thbii),pu,cm,opt.ifinv);                            
      cvcid =covthosl(mean(thci),pm,cm,opt.ifinv);                              
      cvciid =covthosl(mean(thcii),pm,cm,opt.ifinv);                            
      cvdid =covthosl(mean(thdi),pa,cm,opt.ifinv);                              
      cvdiid =covthosl(mean(thdii),pa,cm,opt.ifinv);                            
      cveid =covthosl(mean(thei),pu,cm,opt.ifinv);                              
      cveiid =covthosl(mean(theii),pu,cm,opt.ifinv);                            
      cvfid =covthosl(mean(thfi),pu,cm,opt.ifinv);                              
      cvfiid =covthosl(mean(thfii),pu,cm,opt.ifinv);                            
      cvgid =covthosl(mean(thgi),pu,cm,opt.ifinv);                              
      cvgiid =covthosl(mean(thgii),pu,cm,opt.ifinv);                            
      cvhid =covthosl(mean(thhi),pu,cm,opt.ifinv);                              
      cvhiid =covthosl(mean(thhii),pu,cm,opt.ifinv);                            
      cviid =covthosl(mean(thii),psm,cm,opt.ifinv);                             
      cviiid =covthosl(mean(thiii),psm,cm,opt.ifinv);                           
      cvjid =covthosl(mean(thji),psa,cm,opt.ifinv);                             
      cvjiid =covthosl(mean(thjii),psa,cm,opt.ifinv);  
    end
    % Choose a plot style: style 0 combines the histograms by set of Matern
    % parameters into two rows, style 1 displays many histograms by set of
    % experiment with many rows.
    pltstyle=2;
    if pltstyle==1; numexp=6; else; numexp=10; end
    % Plot the results, first a large row of histograms with two smaller rows of
    % images of the experiments below (pltstyle 0), or as five columns (1)
    thlab={'\sigma^2','\nu','\rho'};
    lets={'a','b','c','d','e','f','g','h','i','j'};
    labs={'Full (un)','Mixed (un)','Mask (ma)','Anti (an)',...
         'Mask (un)','Anti (un)','SmMask (un)','SmAnti(un)',...
         'SmMask (ma)','SmAnti (ma)'};
    labfs={'Full','Mixed','Mask','Anti',...
         'Mask','Anti','Sm.Mask','Sm.Anti','Sm.Mask','Sm.Anti'};
    labts={'(Unit)','(Unit)','(Mask)','(Anti)',...
         '(Unit)','(Unit)','(Unit)','(Unit)','(Sm.Mask)','(Sm.Anti)'};
    lnami ={'(a) I','(b) I','(c) I','(d) I','(e) I','(f) I','(g) I','(h) I','(i) I','(j) I'};
    lnamii={'(a)II','(b)II','(c)II','(d)II','(e)II','(f)II','(g)II','(h)II','(i)II','(j)II'};
    % Plot style 2 order of appearance:
    % estimates
    ps2tind=[1 2 5 6 7 8 3 4 9 10];
    % covariance calculations
    ps2cind=ps2tind;
    % fields
    ps2find=ps2tind;
    % tapers
    ps2tind=ps2tind;%[1 2 5 6 7 8 3 4 9 10];
    if pltstyle==0
      clrs1=spring(numexp);
      clrs2=winter(numexp);
    else
      clrc1=kelicol;
      clrc2=parula;
      clrs1=repmat(clrc1(floor(size(clrc1,1)/2),:),numexp,1);
      clrs2=repmat(clrc2(floor(size(clrc2,1)/2),:),numexp,1);
    end
    numbin=ceil((1+log2(2*numreals*numexp))*4);
    this={thai thbi thci thdi thei thfi thgi thhi thii thji};
    thiis={thaii thbii thcii thdii theii thfii thgii thhii thiii thjii};
    cvis={cvaie cvbie cvcie cvdie cveie cvfie cvgie cvhie cviie cvjie};
    cviis={cvaiie cvbiie cvciie cvdiie cveiie cvfiie cvgiie cvhiie cviiie cvjiie};
    cv1d ={cvaid  cvbid  cvcid  cvdid  cveid  cvfid  cvgid  cvhid cviid cvjid};
    cv2d={cvaiid cvbiid cvciid cvdiid cveiid cvfiid cvgiid cvhiid cviiid cvjiid};
    clf;
    if pltstyle==0
      [ah,ha,H]=krijetem(subnum(4,numexp));
      % The histograms
      for ind=1:np
        axes(H(1,ind))
        hold on
        binedg1=minmax([thai(:,ind);  thbi(:,ind);  thci(:,ind);  thdi(:,ind);...
                        thei(:,ind);  thfi(:,ind);  thgi(:,ind);  thhi(:,ind);...
                        thii(:,ind); thji(:,ind)]);
        binedg2=minmax([thaii(:,ind); thbii(:,ind); thcii(:,ind); thdii(:,ind);...
                        theii(:,ind); thfii(:,ind); thgii(:,ind); thhii(:,ind);...
                        thiii(:,ind); thjii(:,ind)]);
        binv1=linspace(binedg1(1),binedg1(2),numbin);
        binv2=linspace(binedg2(1),binedg2(2),numbin);
        nlin1=linspace(binedg1(1)*0.5,binedg1(2)*2,numreals);
        nlin2=linspace(binedg2(1)*0.5,binedg2(2)*2,numreals);
        for jnd=1:numexp
          % Make the histogram (I)
          hi {ind}(jnd)=histogram(this{jnd}(:,ind) ,binv1,'FaceAlpha',0.5,...
                        'FaceColor',clrs1(jnd,:),'DisplayName',lnami{jnd});
          %hi {ind}(jnd)=histogram(this{jnd}(:,ind) ,numbin,'FaceAlpha',0.5,...
          %              'FaceColor',clrs1(jnd,:),'DisplayName',lnami{jnd});
          % Make the pdf curves from the analytical covariances (II)
          % npd=normpdf(nlin2,nanmean(thiis{jnd}(:,ind)),sqrt(cv2d{jnd}(ind,ind)));
          % npdii{ind}(jnd)=plot(nlin2,npd.*numreals,'Color',clrs2(jnd,:));%,'Alpha',0.5);
          % npd=normpdf(nlin1,nanmean(this {jnd}(:,ind)),sqrt(cv1d {jnd}(ind,ind)));
          % npdi {ind}(jnd)=plot(nlin1,npd.*numreals,'Color',clrs1(jnd,:));%,'Alpha',0.5);
        end
        box on
        longticks
        if ind==1
          lg1=legend();
          ylabel('Count')
        end
        xlm1{ind}=xlim;
        xli1{ind}=xline(th1(ind),'HandleVisibility','off');
        xlii1{ind}=xline(th2(ind),'HandleVisibility','off');

        axes(H(2,ind))
        hold on
        % Keeps the legend organized
        for jnd=1:numexp
          % Make the histogram (II)
          hii{ind}(jnd)=histogram(thiis{jnd}(:,ind),binv2,'FaceAlpha',0.5,...
                        'FaceColor',clrs2(jnd,:),'DisplayName',lnamii{jnd});
          % Make the pdf curves from the empirical covariances (II)
          % npe=normpdf(nlin2,nanmean(thiis{jnd}(:,ind)),sqrt(cviis{jnd}(ind,ind)));
          % npeii{ind}(jnd)=plot(nlin2,npe.*numreals,'Color',clrs2(jnd,:));%,'Alpha',0.5);
          % npe=normpdf(nlin1,nanmean(this {jnd}(:,ind)),sqrt(cvis {jnd}(ind,ind)));
          % npei {ind}(jnd)=plot(nlin1,npe.*numreals,'Color',clrs1(jnd,:));%,'Alpha',0.5);
        end
        box on
        longticks
        xli2{ind}=xline(th1(ind),'HandleVisibility','off');
        xlii2{ind}=xline(th2(ind),'HandleVisibility','off');
        xlm2{ind}=xlim;
        axlm{ind}=[xlm1{ind} xlm2{ind}];
        axlm{ind}=[binedg1 binedg2];
        set(gca,'XTick',[th1(ind) th2(ind) max(axlm{ind})],...
            'XTickLabel',...
            {sprintf('\\theta_I=%0.1f',th1(ind)); sprintf('\\theta_{II}=%0.1f',th2(ind));...
            round(max(axlm{ind}),1)},'XLim',minmax(axlm{ind}))
            % {round(min(axlm{ind}),1); '\theta_I'; '\theta_{II}'; round(max(axlm{ind}),1)},...
            % 'XLim',minmax(axlm{ind}))
        if ind==1
          lg2=legend();
          ylabel('Count')
        end

        % go back and tidy up
        axes(H(1,ind))
        set(gca,'XTick',[th1(ind) th2(ind) max(axlm{ind})],...
            'XTickLabel',...
            {sprintf('\\theta_I=%0.1f',th1(ind)); sprintf('\\theta_{II}=%0.1f',th2(ind));...
            round(max(axlm{ind}),1)},'XLim',minmax(axlm{ind}))
        % set(gca,'XTick',[min(axlm{ind}) th1(ind) th2(ind) max(axlm{ind})],...
            % 'XTickLabel',...
            % {round(min(axlm{ind}),1); '\theta_I'; '\theta_{II}'; round(max(axlm{ind}),1)},...
            % 'XLim',minmax(axlm{ind}))
        tih(ind)=title(sprintf('$%s$',thlab{ind}),'Interpreter','latex');
      end
      % delete the 4th through 8th column of the first row
      delete([ah(4) ah(5) ah(6) ah(7) ah(8)])
      % delete the 2nd row entirely
      %delete([ah(9) ah(10) ah(11) ah(12) ah(13) ah(14) ah(15) ah(16)])
      delete([ah(12) ah(13) ah(14) ah(15) ah(16)])
      % Show the last set of experiment simulations
      Is ={Hai  Hbi  Hci  Hdi  Hci  Hdi  Hgi  Hhi};
      IIs={Haii Hbii Hcii Hdii Hcii Hdii Hgii Hhii};
      for ind=1:numexp
        axes(H(3,ind))
        himg1(ind)=imagesc(v2s(Is{ind},params));
        set(himg1(ind),'AlphaData',~isnan(v2s(IIs{ind},params)))
        if ismember(ind,[3 4 9 10])
          hold on;plot(cr(:,1),cr(:,2),'Color','b','LineWidth',3)
        else
          hold on;
          r(ind)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','b','LineWidth',3);
        end
        axis image
        longticks
        ti3(ind)=title(sprintf('(%s) I: %s',lets{ind},labs{ind}),'Interpreter','latex');
        set(gca,'YTick',[1 params.NyNx(1)],...
                'XTick',[1 params.NyNx(2)],'XTickLabel',['' ''])
        if ind~=1
          set(gca,'YTickLabel',['' ''])
        end

        axes(H(4,ind))
        himg2(ind)=imagesc(v2s(IIs{ind},params));
        set(himg2(ind),'AlphaData',~isnan(v2s(IIs{ind},params)))
        if ismember(ind,[3 4 9 10])
          hold on;plot(cr(:,1),cr(:,2),'Color','b','LineWidth',3)
        else
          hold on;
          r2(ind)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','b','LineWidth',3);
        end
        axis image
        longticks
        ti4(ind)=title(sprintf('(%s)II: %s',lets{ind},labs{ind}),'Interpreter','latex');
        set(gca,'YTick',[1 params.NyNx(1)],'XTick',[1 params.NyNx(2)])
        if ind~=1
          set(gca,'YTickLabel',['' ''])
        end
      end
      % Expand the histograms to fill the top two rows
      shrink([ah(1) ah(2) ah(3) ah(9) ah(10) ah(11)],0.35,1)
      moveh([ah(1) ah(2) ah(3) ah(9) ah(10) ah(11)],0.05)
      movev([ah(1) ah(2) ah(3) ah(9) ah(10) ah(11)],-0.01)
      movev([H(3,:) H(4,:)],-0.05)
      moveh([ah(3) ah(11)],0.4)
      moveh([ah(2) ah(10)],0.2)
      lg1.Position=[0.0338 0.2792 0.0649 0.1645];
      lg2.Position=[0.0354 0.0568 0.0649 0.1645];
      sti(1)=sgtitle(...
        sprintf('$%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s$%i%s%i%s%i%s%i%s%s',...
        '\theta_{I}=[',th1(1),'\,',th1(2),'\,',th1(3),...
        ']\,|\,\theta_{II}=[',th2(1),'\,',th2(2),'\,',th2(3),']\,|\,',...
        params.NyNx(1).*params.dydx(1),' m x ',...
        params.NyNx(2).*params.dydx(2),'m; dy=',params.dydx(1),'m, dx=',...
        params.dydx(2),'m; ',upper(params.mask)),... 
       'Interpreter','latex');
    elseif pltstyle==2
      % Use the kdes from the samplingcases figure
      % greatest spread for all parameters
      maxcv=max(reshape([max(cat(3,cvis{:}),[],3),max(cat(3,cviis{:}),[],3),...
                max(cat(3,cv1d{:}),[],3),max(cat(3,cv2d{:}),[],3)],[3 3 4]),[],3);
      maxsd=sqrt(diag(maxcv));
      datspc=0.1;
      stdOut=5;
      nutickvals=[1/3 1/2 1 3/2 5/2];
      nuticklabs={'$\frac{1}{3}$','$\frac{1}{2}$','$1$','$\frac{3}{2}$','$\frac{     5}{2}$'};
      L=floor(sqrt(sum((params.NyNx.*params.dydx).^2)));
      clrs=["#EDB120","#0072BD"];
      % face and edge alpha
      fa=0.1;
      % font size
      afs=10;
      % power of 2 for kde
      kdep=10;
      % labels
      labs={'$\hat{\sigma}^2$','$\hat{\nu}$','$\hat{\rho}$'};

      % all of the fields
      Is ={Hai  Hbi  Hci  Hdi  Hgi  Hhi  Hci  Hdi Hgi Hhi};                             
      IIs={Haii Hbii Hcii Hdii Hgii Hhii Hcii Hdii Hgii Hhii}; 

      % reorder the estimates and covariances
      this=this(ps2tind);
      thiis=thiis(ps2tind);
      cv1d=cv1d(ps2cind);
      cv2d=cv2d(ps2cind);
      cvis=cvis(ps2cind);
      cviis=cviis(ps2cind);

      % one row per experiment
      numrows=10;

      % prepare common color axes
      imcax=max(abs(minmax([Is{:} IIs{:}])));
      imcax1=max(abs(minmax([Is{:}])));
      imcax2=max(abs(minmax([IIs{:}])));

      % prepare the figure axes
      figure(3)
      clf
      [ah,ha,H]=krijetem(subnum(numrows,5));
      % 10 experiments
      for jnd=1:numrows
        % first current field
        cf1=Is{jnd};
        % first set of estimates
        thhat1=this{jnd};
        % first mean
        mobs1=mean(thhat1);
        % first predicted s.d.
        covsd1=sqrt(diag(cv1d{jnd}));
        % first observed s.d.
        obssd1=sqrt(diag(cvis{jnd}));
        % second current field
        cf2=IIs{jnd};
        % second set of estimates
        thhat2=thiis{jnd};
        % second mean
        mobs2=mean(thhat2);
        % second predicted s.d.
        covsd2=sqrt(diag(cv2d{jnd}));
        % second observed s.d.
        obssd2=sqrt(diag(cviis{jnd}));

        % FIRST REALIZATION
        axes(ah((jnd-1)*5+1))
        hold on
       % himg1(jnd)=imagesc(v2s(cf1,params));
       % colormap(H(jnd,1),clrc1)
       % caxis([-imcax1 imcax1])
       % set(himg1(jnd),'AlphaData',~isnan(v2s(cf1,params)))
        himg1(jnd)=imagefnan([1 1],params.NyNx([2 1]),...
                   v2s(cf1,params),clrc1,[-imcax1 imcax1],[],[],1);
        % visually mark where the taper is: unit boxcar or along the boundary?
        if jnd>6
          hold on;plot(cr(:,1),cr(:,2),'Color','k','LineWidth',1.5)
          % and add this because box on does not work well here
          r(jnd)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','k','LineWidth',0.2);
        else
          hold on;
          r(jnd)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','k','LineWidth',1.5);
        end
        % title
        if jnd==1
          ti1(jnd)=title('I','Interpreter','latex');
        end
        set(gca,'YTick',[1 2*params.NyNx(1)/5 3*params.NyNx(1)/5 params.NyNx(1)],...
                'YTickLabel',[string(1) sprintf('(%s) %s',lets{jnd},labfs{ps2tind(jnd)})...
                              sprintf('%s',labts{ps2tind(jnd)}) string(params.NyNx(1))],...
                'XTick',[1 params.NyNx(2)],'XTickLabelRotation',0)
        % only label the final row
        if jnd~=numrows
          set(gca,'XTickLabel',['' ''])
        end
        axis ij image
        longticks
        box on
        hold off

        % THE KDES
        for knd=1:3
          axes(ah((jnd-1)*5+1+knd))
          hold on
          % plot a thin dotted line at the truth
          plot([th1(knd) th1(knd)],[0 10],':','LineWidth',0.5,...
               'Color',clrs{1},'HandleVisibility','off')
          plot([th2(knd) th2(knd)],[0 10],':','LineWidth',0.5,...
               'Color',clrs{2},'HandleVisibility','off')
          datspc=mean([th1(knd) th2(knd)]);
          kdelim{jnd,knd}=datspc+[-maxsd(knd)*stdOut maxsd(knd)*stdOut];
          xnorm=datspc+linspace(-maxsd(knd)*stdOut,maxsd(knd)*stdOut,1000);
          % first kde
          [a1,bdens1,c1]=kde(thhat1(:,knd),2^kdep,...
            kdelim{jnd,knd}(1),kdelim{jnd,knd}(2));
          t1=bar(c1,bdens1*obssd1(knd),1,'FaceColor',clrs{1},...
            'EdgeColor',clrs{1},'FaceAlpha',fa,'EdgeAlpha',fa);
          pobs1{jnd,knd}=plot(xnorm,obssd1(knd).*...
                              normpdf(xnorm,mobs1(knd),obssd1(knd)),...
                              'HandleVisibility','off');
          pclc1{jnd,knd}=plot(xnorm,covsd1(knd).*...
                               normpdf(xnorm,mobs1(knd),covsd1(knd)),...
                               'HandleVisibility','off');
          set(pclc1{jnd,knd},'LineWidth',0.5,'Color','k')
          set(pobs1{jnd,knd},'LineWidth',1.0,'Color',clrs{1})

          % second kde
          [a2,bdens2,c2]=kde(thhat2(:,knd),2^kdep,...
            kdelim{jnd,knd}(1),kdelim{jnd,knd}(2));
          % use a common scaling with the first set? let's try...
          t2=bar(c2,bdens2*obssd1(knd),1,'FaceColor',clrs{2},...
            'EdgeColor',clrs{2},'FaceAlpha',fa,'EdgeAlpha',fa);
          pobs2{jnd,knd}=plot(xnorm,obssd1(knd).*...
                              normpdf(xnorm,mobs2(knd),obssd2(knd)),...
                              'HandleVisibility','off');
          pclc2{jnd,knd}=plot(xnorm,covsd1(knd).*...
                               normpdf(xnorm,mobs2(knd),covsd2(knd)),...
                               'HandleVisibility','off');
          set(pclc2{jnd,knd},'LineWidth',0.5,'Color','k')
          set(pobs2{jnd,knd},'LineWidth',1.0,'Color',clrs{2})
          %
          maxx(jnd,knd)=max(abs([c1 c2]));
          maxy(jnd,knd)=max(abs([bdens1*obssd1(knd); bdens2*obssd2(knd)]));
          plot([mean([c1 c2])+stdOut.*maxsd(knd).*[-1 1]],[0 0],'k')
          if jnd~=numrows
            set(gca,'XTickLabel',['' ''])
          end
          longticks
          set(gca,'YTickLabel',[])
          box on
          hold off
        end

        % THE SECOND FIELD
        axes(ah(jnd*5))
        hold on
       % himg2(jnd)=imagesc(v2s(cf2,params));
       % colormap(H(jnd,5),clrc2)
       % caxis([-imcax2 imcax2])
       % set(himg2(jnd),'AlphaData',~isnan(v2s(cf2,params)))
        himg2(jnd)=imagefnan([1 1],params.NyNx([2 1]),...
                   v2s(cf2,params),clrc2,[-imcax2 imcax2],[],[],1);
        if jnd>6
          hold on;plot(cr(:,1),cr(:,2),'Color','k','LineWidth',1.5)
          % and add this because box on does not work well here
          r2(jnd)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','k','LineWidth',0.2);
        else
          hold on;
          r2(jnd)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','k','LineWidth',1.5);
        end
        axis ij image
        if jnd==1
          ti5(jnd)=title('II','Interpreter','latex');
        end
        set(gca,'YTick',[1 params.NyNx(1)],'YAxisLocation','right',...
                'XTick',[1 params.NyNx(2)],'XTickLabelRotation',0)
        % only label the bottom row
        if jnd~=numrows
          set(gca,'XTickLabel',['' ''])
        end
        longticks
        box on
        hold off
        %
        % Now assign a common axis range for all of the kde subplots
        for knd=1:3
          ylspc=1.1;
          datspc=mean([th1(knd) th2(knd)]);
          set(H(:,knd+1),'XLim',[th1(knd)-maxsd(knd)*1 th2(knd)+maxsd(knd)*2],...
              'YLim',max(maxy(:,knd)).*[0 ylspc]);
          sth1=[{sprintf('%0.1f',th1(1))} {sprintf('%0.1f',th1(2))}...
                {sprintf('%0.1fe3     ',th1(3)/1e3)}];
          sth2=[{sprintf('%0.1f',th2(1))} {sprintf('%0.1f',th2(2))}...
                {sprintf('     %0.1fe3',th2(3)/1e3)}];
          if jnd==numrows
            xtcklb={sprintf('%s',sth1{knd});...
                    sprintf('%s',sth2{knd})};
            xlb=sprintf('$%s$',thlab{knd});
          else
            xtcklb={'' '' ''};
            xlb='';
          end
          axes(ah((jnd-1)*5+1+knd))
          xlabel(xlb,'Interpreter','latex')
          set(gca,'XTick',[th1(knd) th2(knd)],...
              'XTickLabel',xtcklb,'XTickLabelRotation',0)
          box on
          if jnd==1
            xtcklb1={sprintf('%s',sth1{knd});...
                     sprintf('%s',sth2{knd})};
            axes(H(1,1+knd))
            set(gca,'XTick',[th1(knd) th2(knd)],'XTickLabel',xtcklb1,...
                'XTickLabelRotation',0,'XAxisLocation','top')
          %  [axx(knd),xl(knd),yl(knd)]=...
          %      xtraxis(H(jnd,knd+1),...
          %              [th1(knd) th2(knd)],xtcklb1);
          %  longticks
            tih(knd)=title(sprintf('$%s$',thlab{knd}),'Interpreter','latex');
            box on
          end
        end
      end
      % Add a single legend
      axes(H(1,4))
      [lg,lgi]=legend({'$\theta_I$','$\theta_{II}$'},'Interpreter','latex');
      lg.AutoUpdate='off';
      
      % Copy the inset for the mixed fields:
      % I to II
      tmpfig=figure();
      Is{3}(Is{3}==0)=deal(NaN);
      tmpah1=subplot(1,1,1);
      tmpim1=imagesc(v2s(Is{3},params));
      axis image ij
      colormap(tmpah1,clrc1)
      caxis([-imcax1 imcax1])
      set(tmpim1,'AlphaData',~isnan(v2s(Is{3},params)),'AlphaDataMapping','none')
      cah1=copyobj(tmpah1,tmpfig);
      set(cah1,'Parent',figure(3),'Position',get(H(2,5),'position'),...
               'Color','none','YTick','','XTick','')
      close(tmpfig)
      % II to I
      tmpfig=figure();
      IIs{3}(IIs{3}==0)=deal(NaN);
      tmpah2=subplot(1,1,1);
      tmpim2=imagesc(v2s(IIs{3},params));
      axis image ij
      colormap(tmpah2,clrc2)
      caxis([-imcax2 imcax2])
      set(tmpim2,'AlphaData',~isnan(v2s(IIs{3},params)),'AlphaDataMapping','none')
      cah2=copyobj(tmpah2,tmpfig);
      set(cah2,'Parent',figure(3),'Position',get(H(2,1),'position'),...
               'Color','none','YTick','','XTick','')
      close(tmpfig)
      %
      figure(3)
      sti(1)=sgtitle(...
        sprintf('$%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s$%i%s%i%s%i%s%i%s%s',...
        '\theta_{I}=[',th1(1),'\;',th1(2),'\;',th1(3),...
        ']\;|\;\theta_{II}=[',th2(1),'\;',th2(2),'\;',th2(3),']\;|\;',...
        params.NyNx(1).*params.dydx(1)./1e3,'km x ',...
        params.NyNx(2).*params.dydx(2)./1e3,'km; dy=',params.dydx(1)/1e3,'km, dx=',...
        params.dydx(2)/1e3,'km; ',upper(params.mask)),...
       'Interpreter','latex');

      shrink(H(:,[2:4]),0.75,1)
      moveh(H(:,2),-0.02)
      moveh(H(:,4),0.02)
      movev(H,-0.03)
      movev(cah1,-0.03)
      movev(cah2,-0.03)
      keyboard

      axes(H(1,1))
      [cb(1),xcb(1)]=addcb('vert',[],[],'kelicol',[],0);
      %xcb(1).String='field [units]';
      xcb(1).String='';
      cb(1).TickDir='out';
      cb(1).TickLength(1)=0.01;
      cb(1).YAxisLocation='right';
      cb(1).YTick=linspace(cb(1).YTick(1),cb(1).YTick(end),5);
      hvr=minmax(Hai);
      cb(1).YTickLabel=round(linspace(hvr(1),hvr(2),5),1);
      cb(1).Position(1)=0.225;

      axes(H(1,5))
      [cb(2),xcb(2)]=addcb('vert',[],[],'parula',[],0);
      %xcb(2).String='field [units]';
      xcb(2).String='';
      cb(2).TickDir='out';
      cb(2).TickLength(1)=0.01;
      cb(2).YAxisLocation='left';
      cb(2).YTick=linspace(cb(2).YTick(1),cb(2).YTick(end),5);
      hvr=minmax(Haii);
      cb(2).YTickLabel=round(linspace(hvr(1),hvr(2),5),1);
      cb(2).Position(1)=0.812;

      % For the table
       tnd=2;
       myths=this;
       mycvd=cv1d;
       mycvs=cvis;
       myth0s=th1;
       mobs=mean(myths{tnd});
       sclc=sqrt(diag(mycvd{tnd}));
       sobs=sqrt(diag(mycvs{tnd}));
       mobs(1), sclc(1)/myth0s(1)*100, sobs(1)/myth0s(1)*100
       mobs(2), sclc(2)/myth0s(2)*100, sobs(2)/myth0s(2)*100
       sprintf('%0.2f',mobs(3))
       pi*mobs(3)/L, sclc(3)/myth0s(3)*100, sobs(3)/myth0s(3)*100
       mycvd{tnd}./[diag(sqrt(mycvd{tnd}))*diag(sqrt(mycvd{tnd}))']*100
       mycvs{tnd}./[diag(sqrt(mycvs{tnd}))*diag(sqrt(mycvs{tnd}))']*100
      %
      % Save the figure
      saveas(gcf,...
        sprintf('maskit_demo5_%s_%s.eps',params.mask,date),'epsc')
    elseif pltstyle==1
      % We will make this figure twice -- the first time for the masked fields,
      % the second time for the anti-masked fields
      for fnd=1:2
      figure(fnd)
      if fnd==1
        % the masks
        Inds=[1 2 3 5 7 9];
      else
        % the anti-masks
        Inds=[1 2 4 6 8 10];
      end
      [ah,ha,H]=krijetem(subnum(numexp,5));
      % Realizations of the experiments as two columns (number of models)
      % Show the last set of experiment simulations
      allIs ={Hai  Hbi  Hci  Hdi  Hci  Hdi  Hgi  Hhi Hgi Hhi};
      allIIs={Haii Hbii Hcii Hdii Hcii Hdii Hgii Hhii Hgii Hhii};
      Is ={allIs{Inds}};
      IIs={allIIs{Inds}};
      let={lets{Inds}};
      labf={labfs{Inds}};
      labt={labts{Inds}};

      imcax=max(abs(minmax([Is{:} IIs{:}])));
      imcax1=max(abs(minmax([Is{:}])));
      imcax2=max(abs(minmax([IIs{:}])));
      for ind=1:numexp
        % Experiment realization for I
        axes(H(ind,1))
        himg1(ind)=imagesc(v2s(Is{ind},params));
        colormap(H(ind,1),clrc1)
        caxis([-imcax1 imcax1])
        set(himg1(ind),'AlphaData',~isnan(v2s(Is{ind},params)))
        if ismember(ind,[3 6])
          hold on;plot(cr(:,1),cr(:,2),'Color','b','LineWidth',3)
        else
          hold on;
          r(ind)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','b','LineWidth',3);
        end
        axis image
        longticks
        if ind==1
          ti3(ind)=title('I','Interpreter','latex');
        end
        set(gca,'YTick',[1 2*params.NyNx(1)/5 3*params.NyNx(1)/5 params.NyNx(1)],...
                'YTickLabel',[string(1) sprintf('(%s) %s',let{ind},labf{ind})...
                              sprintf('%s',labt{ind}) string(params.NyNx(1))],...
                'XTick',[1 params.NyNx(2)],'XTickLabelRotation',0)
        if ind~=numexp
          set(gca,'XTickLabel',['' ''])
        end

        % Experiment realization for II
        axes(H(ind,2))
        himg2(ind)=imagesc(v2s(IIs{ind},params));
        colormap(H(ind,2),clrc2)
        caxis([-imcax2 imcax2])
        set(himg2(ind),'AlphaData',~isnan(v2s(IIs{ind},params)))
        if ismember(ind,[3 6])
          hold on;plot(cr(:,1),cr(:,2),'Color','b','LineWidth',3)
        else
          hold on;
          r2(ind)=rectangle('Position',[1 1 params.NyNx(2)-0.5 params.NyNx(1)-0.5],...
                            'EdgeColor','b','LineWidth',3);
        end
        axis image
        longticks
        if ind==1
          ti4(ind)=title('II','Interpreter','latex');
        end
        set(gca,'YTick',[1 params.NyNx(1)],'YTickLabel',['' ''],...
                'XTick',[1 params.NyNx(2)],'XTickLabelRotation',0)
        if ind~=numexp
          set(gca,'XTickLabel',['' ''])
        end
      end

      thisf ={this{Inds}};
      thiisf={thiis{Inds}};
      % The histograms as 3 columns (number of Matern parameters)
      for ind=1:np
        binedg=[cat(1,thisf{:}); cat(1,thiisf{:})];
        binedg=minmax(binedg(:,ind)).*[0.5 1.2];
        binv=linspace(binedg(1),binedg(2),numbin);
        axlm{ind}=[];
        % Generate normally distributed samples for the parameter and its
        % standard deviation calculated from the analytical covariance for the
        % parameter value and its taper
        dmy1a=random('normal',mean(thai(:,ind)), sqrt(cvaid(ind,ind)) ,numreals,1);
        dmy1b=random('normal',mean(thbi(:,ind)), sqrt(cvbid(ind,ind)) ,numreals,1);
        dmy1c=random('normal',mean(thci(:,ind)), sqrt(cvcid(ind,ind)) ,numreals,1);
        dmy1d=random('normal',mean(thdi(:,ind)), sqrt(cvdid(ind,ind)) ,numreals,1);
        dmy1e=random('normal',mean(thei(:,ind)), sqrt(cveid(ind,ind)) ,numreals,1);
        dmy1f=random('normal',mean(thfi(:,ind)), sqrt(cvfid(ind,ind)) ,numreals,1);
        dmy1g=random('normal',mean(thgi(:,ind)), sqrt(cvgid(ind,ind)) ,numreals,1);
        dmy1h=random('normal',mean(thhi(:,ind)), sqrt(cvhid(ind,ind)) ,numreals,1);
        dmy1i=random('normal',mean(thii(:,ind)), sqrt(cviid(ind,ind)) ,numreals,1);
        dmy1j=random('normal',mean(thji(:,ind)), sqrt(cvjid(ind,ind)) ,numreals,1);
        dmy2a=random('normal',mean(thaii(:,ind)),sqrt(cvaiid(ind,ind)),numreals,1);
        dmy2b=random('normal',mean(thbii(:,ind)),sqrt(cvbiid(ind,ind)),numreals,1);
        dmy2c=random('normal',mean(thcii(:,ind)),sqrt(cvciid(ind,ind)),numreals,1);
        dmy2d=random('normal',mean(thdii(:,ind)),sqrt(cvdiid(ind,ind)),numreals,1);
        dmy2e=random('normal',mean(theii(:,ind)),sqrt(cveiid(ind,ind)),numreals,1);
        dmy2f=random('normal',mean(thfii(:,ind)),sqrt(cvfiid(ind,ind)),numreals,1);
        dmy2g=random('normal',mean(thgii(:,ind)),sqrt(cvgiid(ind,ind)),numreals,1);
        dmy2h=random('normal',mean(thhii(:,ind)),sqrt(cvhiid(ind,ind)),numreals,1);
        dmy2i=random('normal',mean(thiii(:,ind)),sqrt(cviiid(ind,ind)),numreals,1);
        dmy2j=random('normal',mean(thjii(:,ind)),sqrt(cvjiid(ind,ind)),numreals,1);

        dmy1s={dmy1a dmy1b dmy1c dmy1d dmy1e dmy1f dmy1g dmy1h dmy1i dmy1j};
        dmy2s={dmy2a dmy2b dmy2c dmy2d dmy2e dmy2f dmy2g dmy2h dmy2i dmy2j};
        dmy1={dmy1s{Inds}};
        dmy2={dmy2s{Inds}};
        for jnd=1:numexp
          axes(H(jnd,np+ind-1))
          hold on
          % Make a histogram for each experiment (I), controlling for the number
          % and placement of bins
          hi {ind}(jnd)=histogram(thisf{jnd}(:,ind) ,binv,'FaceAlpha',0.5,...
                        'FaceColor',clrs1(jnd,:));
          % Make another histogram (II)
          hii{ind}(jnd)=histogram(thiisf{jnd}(:,ind),binv,'FaceAlpha',0.5,...
                        'FaceColor',clrs2(jnd,:));
          % We will place a single legend in the top right subplot
          if jnd==1&ind==3
            % Include a marker for predicted distributions?
            lg=legend({'$\hat{\theta}_I$','$\hat{\theta}_{II}$'});
            lg.Interpreter='latex';
            lg.AutoUpdate='off';
          end

          % Trick histfit into just plotting a curve for the empirical
          % distribution
          hdta=(hi{ind}(jnd).BinEdges(1:end-1)+hi{ind}(jnd).BinEdges(2:end))/2;
          hdtav=repelem(hdta,hi{ind}(jnd).Values);
          fhi{ind}(jnd,:)=histfit(hdtav,sum(hi{ind}(jnd).Values>0));
          set(fhi{ind}(jnd,2),'Color','#EDB120','LineWidth',1.5)
          delete(fhi{ind}(jnd,1))
          hdta=(hii{ind}(jnd).BinEdges(1:end-1)+hii{ind}(jnd).BinEdges(2:end))/2;
          hdtav=repelem(hdta,hii{ind}(jnd).Values);
          fhii{ind}(jnd,:)=histfit(hdtav,sum(get(hii{ind}(jnd),'Values')>0));
          set(fhii{ind}(jnd,2),'Color','#0072BD','LineWidth',1.5)
          delete(fhii{ind}(jnd,1))

          % Make curves from the analytical covariances (II) by sampling
          % the normal distribution and using the histfit trick from above
          % npd=pdf('normal',th1(ind),sqrt(cv1d{jnd}(ind,ind)));
          %nbn=sum(get(hi{ind}(jnd),'Values')>0);
          %mn=th1(ind); st=sqrt(cv1d{jnd}(ind,ind));
          %dmy = max(nbn)*exp(-(((nlin-mn).^2)/(2*st.^2)));
          %plot(nlin,dmy,'k','LineWidth',1)
          dmyh=histogram(dmy1{jnd},binv);
          fdta=(dmyh.BinEdges(1:end-1)+dmyh.BinEdges(2:end))/2;
          fdtav=repelem(fdta,dmyh.Values);
          cvhi{ind}(jnd,:)=histfit(fdtav,sum(dmyh.Values>0));
          delete(cvhi{ind}(jnd,1))
          delete(dmyh)
          set(cvhi{ind}(jnd,2),'Color','k','LineWidth',0.5,'LineStyle','-')
          
          % And the second one, too
          % npd=pdf('normal',th2(ind),sqrt(cv2d{jnd}(ind,ind)));
          %nbn=sum(get(hi{ind}(jnd),'Values')>0);
          %mn=th2(ind); st=sqrt(cv2d{jnd}(ind,ind));
          %dmy = max(nbn)*exp(-(((nlin-mn).^2)/(2*st.^2)));
          %% dmy = max(numbin)*exp(-(((nlin-npd.mean).^2)/(2*npd.sigma.^2)));
          %plot(nlin,dmy,'k','LineWidth',1)
          dmyh=histogram(dmy2{jnd},binv);
          fdta=(dmyh.BinEdges(1:end-1)+dmyh.BinEdges(2:end))/2;
          fdtav=repelem(fdta,dmyh.Values);
          cvhii{ind}(jnd,:)=histfit(fdtav,sum(dmyh.Values>0));
          delete(cvhii{ind}(jnd,1))
          delete(dmyh)
          set(cvhii{ind}(jnd,2),'Color','k','LineWidth',0.5,'LineStyle','-')

          xli2{jnd}=xline(th1(ind),'HandleVisibility','off');
          xlii2{jnd}=xline(th2(ind),'HandleVisibility','off');
          if ind==1
            ylabel('Count')
          else
            set(gca,'YTickLabel','')
          end
          box on
          longticks
          xlm=xlim;
          axlm{ind}=[axlm{ind} {xlim}];
        end
        axlm{ind}=minmax([axlm{ind}{:}]);
        % We don't want to see too far beyond where the data appears
        axlm{ind}(axlm{ind}<0)=0;
        axlm{ind}(axlm{ind}<0.6*binedg(1))=0.6*binedg(1);
        axlm{ind}(axlm{ind}>1.05*binedg(2))=1.05*binedg(2);
        % Go back through and adjust xlim
        for jnd=1:numexp
          axes(H(jnd,ind+np-1))
          set(gca,'XLim',axlm{ind})
          sth1=[{sprintf('%0.1f',th1(1))} {sprintf('%0.1f',th1(2))}...
                {sprintf('%0.1fe3  ',th1(3)/1e3)}]; 
          sth2=[{sprintf('%0.1f',th2(1))} {sprintf('%0.1f',th2(2))}...
                {sprintf('  %0.1fe3',th2(3)/1e3)}]; 
          %xtcklb={round(min(axlm{ind}),1); sprintf('\\theta_I=%s',sth1{ind});...
          if ind==3
            ltcklb=sprintf('%0.1fe3',round(max(axlm{ind})/1e3,1));
          else
            ltcklb=round(max(axlm{ind}),1);
          end
          xtcklb={sprintf('%s',sth1{ind});...
            sprintf('%s',sth2{ind});...
            ltcklb};
          % xtcklb={sprintf('\\theta_I=%s',sth1{ind});...
          %   sprintf('\\theta_{II}=%s',sth2{ind});...
          %   ltcklb};
          if jnd==numexp
            mxtcklb=xtcklb;
          else
            %mxtcklb={'' '' '' ''};
            mxtcklb={'' '' ''};
          end
          %set(gca,'XTick',[min(axlm{ind}) th1(ind) th2(ind) max(axlm{ind})],...
          set(gca,'XTick',[th1(ind) th2(ind) max(axlm{ind})],...
              'XTickLabel',mxtcklb,'YLim',[0 numreals],'XTickLabelRotation',0)
              %'XLim',axlm{ind},'YLim',[0 numreals])
          if jnd==1
            [axx(ind),xl(ind),yl(ind)]=...
              xtraxis(H(jnd,np+ind-1),...
                      [th1(ind) th2(ind) max(axlm{ind})],xtcklb);
                      %[min(axlm{ind}) th1(ind) th2(ind) max(axlm{ind})],xtcklb);
             longticks
             tih(ind)=title(sprintf('$%s$',thlab{ind}),'Interpreter','latex');
          end
        end
      end
      sti(1)=sgtitle(...
        sprintf('$%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s$%i%s%i%s%i%s%i%s%s',...
        '\theta_{I}=[',th1(1),'\;',th1(2),'\;',th1(3),...
        ']\;|\;\theta_{II}=[',th2(1),'\;',th2(2),'\;',th2(3),']\;|\;',...
        params.NyNx(1).*params.dydx(1)./1e3,'km x ',...
        params.NyNx(2).*params.dydx(2)./1e3,'km; dy=',params.dydx(1)/1e3,'km, dx=',...
        params.dydx(2)/1e3,'km; ',upper(params.mask)),... 
       'Interpreter','latex');
      moveh(H(:,[1 2]),-0.1)
      moveh(H(:,2),-0.04)
      shrink([H(:,[1:2])],0.62, 0.82)
      shrink([H(:,[3:5]); axx(:)'],0.62,0.82)
      moveh([H(:,3); axx(1)],-0.11)
      moveh([H(:,4); axx(2)],-0.05)
      moveh([H(:,5); axx(3)],0.02)
      movev(H,-0.04); movev([axx(:)],-0.04)
      moveh(H,0.03); moveh([axx(:)],0.03)
      moveh([H(:,[3:5]); axx(:)'],-0.03)
      moveh(H(:,2),-0.05)
      movev([tih(:)],-5)
    end
    keyboard
    % Save the figure
    saveas(gcf,...
      sprintf('maskit_demo5_%s_%s_%i.eps',params.mask,date,ind),'epsc')
    end
    % To get around the output wanting more
    w=NaN; vw=NaN; cr=NaN;
end

% Variable output
varns={v,I,w,vw,scl,cr};
varargout=varns(1:nargout);

function plotit(v,p,cr,th)
% Should use th(1) instead of halverange, shouldn't I...
imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],halverange(v,80,NaN)); 
if strcmp(p.mask,'france')
  axis ij image 
else
  axis xy image
end
hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
% We might be bringing in a mask that has been shifted off of the main grid and
% cr will plot beyond the boundaries of the field. Set axes limits for the
% visualization
ylim([1 p.NyNx(1)]); xlim([1 p.NyNx(2)])
if ~isempty(th)
    if isstr(th)
        t=title(th);
    else
        t=title(sprintf('%i x %i | %s = [%g %g %gx]',p.NyNx(1),p.NyNx(2),'\theta',...
                        th./[1 1 sqrt(prod(p.dydx))]));
    end
    movev(t,-p.NyNx(1)/20)
end
