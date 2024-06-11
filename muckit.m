function varargout=muckit(v,p,scl,opt)
% [v,I,w,scl,cr]=MUCKIT(v,p,scl,opt)
%
% Makes randomly mucked up speckled indicatrix tapers.
%
% INPUT:
%
% v       A vector that is the unwrapping of a matrix
% p       A single parameter structure with AT LEAST this:
%           NyNx  number of samples in the y and x directions
%           mask  an index matrix with 1 for yes and 0 for no, OR
%                 a muck type name, e.g. 'random' [default]
% scl     A scale [0 to 1] expressing the non-missing data fraction
% opt     A field structure for optimization flags:
%           algo  the optimization algorithm
%           ifinv  ordered inversion flags for [s2 nu rho], e.g. [1 0 1] 
%
% OUTPUT:
%
% v       The first input with the samples in I retained, rest set to NaN
% I       The mask, unrolled into a vector, note, this is an "anti-mask"
% w       The antimucked-up output matrix, unrolled into a vector
% scl     A scale [0 to 1] expressing the non-missing data fraction
% cr      The colum,row index coordinates of the bounding curve
%
% EXAMPLE:
%
% [~,I]=muckit('demo1','random',0.8); % To get a quick muck to look at
% muckit('demo1') % A uniformly random masking a single field
% muckit('demo1',[],rand) % for a different percentage of survivors
% muckit('demo2'); % To do a serious of mock muck inversions
% maskit('demo3',[],rand) % Illustrating the nomenclature
%
% Last modified by owalbert-princeton.edu, 06/11/2024
% Last modified by fjsimons-at-alum.mit.edu, 06/11/2024

% The default is the demo, for once
defval('v','demo1')
defstruct('opt',{'algo','ifinv'},{'unc',[1 1 1]});

if ~isstr(v)
    % Get the muck by name? Let's play "mikado"?
    if isstr(p.mask)

        % Scale indicates survivorship (in symmetry with MASKIT)
        % Default scale so bounding-box-area-proportionality is scl
        defval('scl',0.9)

	if strcmp(p.mask,'random')
	  I=zeros(p.NyNx);
	  % Now knock out the missing data
	  % This produces way too much overlap
	  % round(scale(rand(round(scl*prod(p.NyNx)),1),[1 prod(p.NyNx)]))
	  % This is exact
	  S=shuffle(1:prod(p.NyNx));
	  I(S(1:round(scl*prod(p.NyNx))))=1;
	  % Produce the bounding curve (later)
	  cr=NaN;
	end
    else
      % You already have it
      cr=NaN;
      I=p.mask;
    end
    w=v;
    % Apply the muck to the field
    v(~I(:))=NaN;
    % Apply the antimuck to the field
    w(~~I(:))=NaN;
    % Note that you must not feed fields with NaN into MLEOSL - the taper I will
    % take of the selectivity
elseif strcmp(v,'demo1')
    % Capture the second input
    defval('p','random')
    defp=p; clear p
    % Now proceed with a fresh copy
    p.mask=defp;
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc; 
    [Hx,th,p]=simulosl([],p,1);
    % Capture the third input, default is none, which defaults inside
    defval('scl',[])
    [Hm,I,~,scl,cr]=muckit(Hx,p,scl);
    clf
    ah(1)=subplot(121); plotit(Hx,p,cr,th)
    ah(2)=subplot(122); plotit(Hm,p,cr,[])
    t=title(sprintf('%i %% surviving',round(100*scl)));
    movev(t,-p.NyNx(1)/20)
    % Create right type of output even though there is non
    v=Hm; w=NaN;
elseif strcmp(v,'demo2')
    % Capture the second input
    defval('p','random')
    defp=p; clear p
    % Capture the fourth input
    defstruct('opt',{'algo','ifinv'},{[],[]});

    % Change the below if you like
    opt.ifinv=[1 1 1];

    % Now proceed with a fresh copy
    p.mask=defp;
    % Simulate using circulant embedding, no taper
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc;
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);
    % Something larger without overdoing it, check weirdness
    % If you want to generate the same hash you need to keep the same dimensions
    p.NyNx=[196 238];

    % Do all the tests or not
    xver=0;
    % Decide on the length of the pause between displays, if any, 
    if xver==1
        pz=2;
    else
        pz=0.1;
    end

    % Define the parameters for the field
    th=[1.50 0.5 30000];

    % Needing to be very explicit in order for parallel computing to work with in the loop
    p.dydx=[10000 10000];
    p.taper=0; p.nugget=0;
    
    % Number of processors, must agree with your machine
    NumWorkers=8;
    % Number of identical experiments to run experiments
    N=1*NumWorkers;

    % Save some output for later? Make new filename hash from all relevant input
    fname=hash([struct2array(orderfields(p)) struct2array(orderfields(opt)) th     N],'SHA-1');
    
    % You need to have an environmental variable file structure set up
    fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s_%s.mat','MUCKIT',fname));

    % Run the experiment anew
    if ~exist(fnams,'file') % || 1==1
        % Initialize the pool of workers
        if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end

        % Run the experiment!
        parfor index=1:N
            clc; disp(sprintf('\n Simulating the field \n'))

            % Simulate the field
            Hx=simulosl(th,p,1);

            % How much should the surviving area occupy?
            % Don't SET this here or PARFOR won't work
            scl=[];

            % Now do the masking and the merging
            [Hm,I,Ham,scl,cr]=muckit(Hx,p,scl);

            % So Ham is the anti-mucked-up field
            v=Hm; w=Ham;

            % As appropriate you'll force the use of BLUROSY in MATERNOSP in LOGLIOSL in MLEOSL
            % This is confusing PARFOR but works with FOR
            % p.blurs=-1;
            % Make a close initial guess?
            thini(index,:)=th+(-1).^randi(2,[1 3]).*th/1000;

            % Don't perturb the thini for the one you don't want to move; could
            % have also have stuck in th          into the aguess slot of
            % MLEOSL, whence they would be appropriately (or not) perturbed
            if opt.ifinv==[1 0 1]
                % This didn't work
                %thini(index,2)=th(2);
                % But this did
                thini(index,:)=th+(-1).^randi(2,[1 3]).*[th(1) 0 th(3)]/1000;
            end

            % Recover the parameters of the full original field without any mucking
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p1=p; p1.taper=0; p1=rmfield(p1,'mask');
            pause(pz); clc; disp(sprintf('\n Estimating first whole field \n'))
            [thhat3(index,:),~,~,scl3(index,:)]=mleosl(Hx,thini(index,:),p1,opt.algo,[],[],opt.ifinv,xver);
            
            if xver==1
                % Explicitly check the likelihood?
	        disp(sprintf('\n Likelihood at the initial guess %g\n',...
	                     logliosl(knums(p1),thini(index,:)./scl3(index,:),scl3(index,:),p1,...
                                      tospec(Hx(:),p1)/(2*pi))));
	        disp(sprintf('\n Likelihood at the truth    %g\n',...
	                     logliosl(knums(p1),th./scl3(index,:),scl3(index,:),p1,...
                                      tospec(Hx(:),p1)/(2*pi))));
	         % --> inside protect by looking at momx?? That's where it goes wrong
            end

            % Now recover the parameters of the speckled field
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p2=p; p2.taper=I; p2.mask='random';
            pause(pz); clc; disp(sprintf('\n Estimating first speckled field \n'))
            [thhat1(index,:),~,~,scl1(index,:)]=mleosl(Hx,thini(index,:),p2,opt.algo,[],[],opt.ifinv,xver);

            if xver==1
	        % Explicitly check the likelihood?
	        disp(sprintf('\n Likelihood at the initial guess %g\n',...
	                     logliosl(knums(p2),thini(index,:)./scl1(index,:),scl1(index,:),p2,...
	                              tospec(p2.taper(:).*Hx(:),p2)/(2*pi)/sqrt(sum(p2.taper(:).^2))*sqrt(prod(p2.NyNx)),1)));
	        disp(sprintf('\n Likelihood at the truth    %g\n',...
	                     logliosl(knums(p2),th./scl1(index,:),scl1(index,:),p2,...
	                              tospec(p2.taper(:).*Hx(:),p2)/(2*pi)/sqrt(sum(p2.taper(:).^2))*sqrt(prod(p2.NyNx)),1)));
            end

            % Now recover the parameters of the complement to the speckled field
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p3=p; p3.taper=~I; p3.mask='random';
            pause(pz); clc; disp(sprintf('\n Estimating anti-speckled field \n'))
            [thhat2(index,:),~,~,scl2(index,:)]=mleosl(Hx,thini(index,:),p3,opt.algo,[],[],opt.ifinv,xver);
            
            % Pause so you can watch live (big pause if xver==0)
            pause(pz)

            % This is unnecessary for PARFOR but works with FOR
            % Reset to the original values
            %p.blurs=Inf; p.taper=0;
        end
        % Save all the output that you'll need later, may add more later
        save(fnams,'th',      'p','thhat1','thhat2','scl1','scl2',...
             'thhat3','scl3',                                'thini'         ,'N')
             
    else
        % Load all the saved output
        load(fnams)
    end

    % And now look at the statistics of the recovery
    str0='%8.2f %5.2f %6.0f';
    str1=sprintf('%s %s\n',str0,str0);
    str2=sprintf('%s\n',str0);
    disp(sprintf('\nWhole | Speckled\n'))
    disp(sprintf(str1,[thhat3.*scl3 thhat1.*scl1]'))

    disp(sprintf('\nWhole | Anti-Speckled\n'))
    disp(sprintf(str2,[thhat3.*scl3 thhat2.*scl2]'))

    % Save the figures bjootifooly
    % Remake one sampled field just to make the visual
    Hx=simulosl(th,p,1);
    [Hm,cr,I,scl]=muckit(Hx,p,[]);

    % Make a visual for good measure
    figure(3)
    clf
    ah(1)=subplot(121); plotit(Hx,p,[],th)
    xlabel('Original field')
    ah(2)=subplot(122); plotit(Hm,p,[],th); 
    t=title(sprintf('%i %% surviving',round(100*scl)));
    movev(t,-p.NyNx(1)/20)

    % Then plot these things using MLEPLOS
    % So this is the region of speckles
    figure(1)
    mleplos(thhat1.*scl1,th,[],[],[],[],[],p,sprintf('%i %% speckle',round(100*scl)),[],opt.ifinv)
    figure(1)
    figdisp(sprintf('%s_2a',pref(sprintf('%s_%s.mat','MUCKIT',fname))),[],[],2)
    clf
    figure(2)
    figdisp(sprintf('%s_2b',pref(sprintf('%s_%s.mat','MUCKIT',fname))),[],[],2)
    clf
    % And this is the antispeckled region
    figure(1)
    mleplos(thhat2.*scl2,th,[],[],[],[],[],p,sprintf('%i %% antispeckle',round(100*(1-scl))),[],opt.ifinv)
    pause(1)
    figure(1)
    figdisp(sprintf('%s_3a',pref(sprintf('%s_%s.mat','MUCKIT',fname))),[],[],2)
    clf
    figure(2)
    figdisp(sprintf('%s_3b',pref(sprintf('%s_%s.mat','MUCKIT',fname))),[],[],2)

    % With PARFOR none of the once-used are available out of the loop but in
    % this demo2 you don't want any output anyway, so put in empties
    [v,I,w,scl,cr]=deal(NaN);
elseif strcmp(v,'demo3')
    th1=[1.15 1.15 23000];
    p.mask='random';
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);
    % Capture the third input, default is none, which defaults inside
    defval('scl',[])

    % Remake one sampled field just to make the visual
    [Hx,~,p,~,~,Sb1]=simulosl(th1,p);
    [v,I,w,scl,cr]=muckit(Hx,p,scl);

    % Make a visual for good measure
    clf
    ah(1)=subplot(221); plotit(v,p,cr,sprintf('%i %% speckle',round(100*scl)))
    ah(3)=subplot(223); plotit(w,p,cr,sprintf('%i %% antispeckle',round(100*(1-scl))))
    ah(4)=subplot(224); plotit(Hx,p,cr,'original')
    movev(ah(4),0.25)
end

% Variable output
varns={v,I,w,scl,cr};
varargout=varns(1:nargout);

function plotit(v,p,cr,th)
% Should use th(1) instead of halverange, shouldn't I...
imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],halverange(v,80,NaN)); axis ij image 

xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
if ~isempty(th)
    if isstr(th)
        t=title(th);
    else
        t=title(sprintf('%i x %i | %s = [%g %g %gx]',p.NyNx(1),p.NyNx(2),'\theta',...
                        th./[1 1 sqrt(prod(p.dydx))]));
    end
    movev(t,-p.NyNx(1)/20)
end
