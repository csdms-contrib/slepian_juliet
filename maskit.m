function varargout=maskit(v,p,scl,w)
% [v,cr,I,w,vw]=MASKIT(v,p,scl,w)
%
% Arbitrary masking of a gridded spatial field
%
% INPUT:
%
% v       A vector that is the unwrapping of a matrix
% p       A single parameter structure with AT LEAST this:
%           NyNx  number of samples in the y and x directions
%           mask an index matrix with a mask, OR:
%                 a region name that will be scaled to fit
% scl     A scale between 0 and 1 with the occupying center fraction
%         of the domain enclosed by the input curve; due to strict
%         inequality at maximum always leaves one pixel rim; and
%         not being used if an actual mask is input instead of a name
% w       A second vector that is the unwrapping of a second matrix
%
% OUTPUT:
%
% v       The masked output matrix, unrolled into a vector
% cr      The colum,row index coordinates of the masking curve
% I       The mask, unrolled into a vector, note, this is an "anti-mask"
% w       The second masked output matrix, unrolled into a vector
% vw      The merged field, second field w sucked into first field v
%
% EXAMPLE:
%
% [~,~,I]=maskit('demo1'); % To get a quick mask to look at
% maskit('demo1') % A geographical region masking a single field
% maskit('demo1','england') % 'amazon', 'orinoco', for geographical variability
% maskit('demo2') % A geographical region merging two fields, with estimation
% maskit('demo2','england') % 'amazon', 'orinoco', for geographical variability
%
% Last modified by fjsimons-at-alum.mit.edu, 06/07/2024

% The default is the demo, for once
defval('v','demo1')

if ~isstr(v)
    % Get the mask by name
    if isstr(p.mask)
        XY=eval(p.mask);

        % Default scale is 80% of the area
        defval('scl',0.8)
        
        % Scale the mask, remember for the curve the Y goes up, for the image, down
        cr(:,1)=          scale(XY(:,1),[1 p.NyNx(2)]+[1 -1]*(1-scl)/2*(p.NyNx(2)-1));
        cr(:,2)=p.NyNx(1)-scale(XY(:,2),[1 p.NyNx(1)]+[1 -1]*(1-scl)/2*(p.NyNx(1)-1))+1;
        
        % The underlying image grid
        [X,Y]=meshgrid(1:p.NyNx(2),1:p.NyNx(1));
        
        % Determine the mask
        I=inpolygon(X,Y,cr(:,1),cr(:,2));
    else
        % You already have it
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
            vw(I(:)) =w(I(:));
            vw(~I(:))=v(~I(:));
        end
    end
    % Apply the mask to the first field
    v(~I(:))=NaN;
elseif strcmp(v,'demo1')
    % Capture the second input
    defval('p','france')
    defp=p; clear p
    % Now proceed with a fresh copy
    p.mask=defp;
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc; 
    [Hx,th,p]=simulosl([],p,1);
    [Hm,cr,I]=maskit(Hx,p);
    clf
    subplot(121); plotit(Hx,p,cr,th)
    subplot(122); plotit(Hm,p,cr,[])
    v=Hm; [w,vw]=deal(NaN);
elseif strcmp(v,'demo2')
    % Capture the second input
    defval('p','france')
    defp=p; clear p
    % Now proceed with a fresh copy
    p.mask=defp;
    % Simulate using circulant embedding, no taper
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc;
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);
    % Something larger without overdoing it, check weirdness
    % Here is one that has failed in the past
    p.NyNx=[201 241];
    
    % Do all the tests or not
    xver=0;
    % Decide on the length of the pause between displays, if any, 
    if xver==1
        pz=2;
    else
        pz=0.1;
    end

    % Define the parameters for the two fields
    th1=[1.50 0.75 30000];
    th2=[1.90 2.25 10000];

    % Needing to be very explicit in order for parallel computing to work with in the loop
    p.dydx=[10000 10000];
    p.taper=0; p.nugget=0;
    
    % Number of processors, must agree with your machine
    NumWorkers=8;
    % Number of identical experimnets to run experiments
    N=10*NumWorkers;

    % Save some output for later? Make new filename hash from all relevant input
    fname=hash([struct2array(orderfields(p)) th1 th2 N],'SHA-1');
    
    % You need to have an environmental variable file structure set up
    fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s_%s.mat','MASKIT',fname));

    % Run the experiment anew
    if ~exist(fnams,'file')
        % Initialize the pool of workers
        if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end

        % Run the experiment!
        parfor index=1:N
            clc; disp(sprintf('\n Simulating all fields \n'))

            % Simulate first field
            Hx=simulosl(th1,p,1);
            % Simulate second field
            Gx=simulosl(th2,p,1);

            % How much should the masked area occupy?
            scl=[];
            
            % Now do the masking and the merging
            [Hm,cr,I,Gm,HG]=maskit(Hx,p,scl,Gx);

            % So HG is the mixed field
            v=Hm; w=Gm; vw=HG;

            % As appropriate you'll force the use of BLUROSY in MATERNOSP in LOGLIOSL in MLEOSL
            % This is confusing PARFOR but works with FOR
            % p.blurs=-1;

            % Make a close initial guess?
            thini1(index,:)=th1+(-1).^randi(2,[1 3]).*th1/1000;
            thini2(index,:)=th2+(-1).^randi(2,[1 3]).*th2/1000;

            % Recover the parameters of the full original fields without any masking
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p1=p; p1.taper=0;
            pause(pz); clc; disp(sprintf('\n Estimating first whole field \n'))
            [thhat4(index,:),~,~,scl4(index,:)]=mleosl(Hx,thini1(index,:),p1,[],[],[],xver);
            pause(pz); clc ; disp(sprintf('\n Estimating second whole field \n'))
            [thhat5(index,:),~,~,scl5(index,:)]=mleosl(Gx,thini2(index,:),p1,[],[],[],xver);

            % Now recover the parameters of the mixed field but only in the region I or ~I
            % Perform the optimization on the complement which should look like the first field
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p2=p; p2.taper=~I;
            pause(pz); clc; disp(sprintf('\n Estimating first partial field \n'))
            [thhat1(index,:),~,~,scl1(index,:)]=mleosl(HG,thini1(index,:),p2,[],[],[],xver);
            
            % Now recover the parameters of the insert which should look like the second field
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p3=p; p3.taper=I;
            pause(pz); clc; disp(sprintf('\n Estimating second partial field \n'))
            [thhat2(index,:),~,~,scl2(index,:)]=mleosl(HG,thini2(index,:),p3,[],[],[],xver);

            % Now recover the parameters of the merged field without knowing of the partition
            % The relabeling is to make PARFOR work, with FOR they could all just be p
            p4=p; p4.taper=0;
            pause(pz); clc; disp(sprintf('\n Estimating whole mixed field \n'))
            [thhat3(index,:),~,~,scl3(index,:)]=mleosl(HG,[],p4,[],[],[],xver);

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
    [~,cr,~,~,HG]=maskit(Hx,p,[],Gx);

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
    ah(1)=subplot(221); plotit(Hx,p,cr,th1)
    ah(3)=subplot(223); plotit(Gx,p,cr,th2)
    ah(4)=subplot(224); plotit(HG,p,cr,[])
    movev(ah(4),0.25)
    figdisp(sprintf('%s_1',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)

    % Then plot the estimates using MLEPLOS, remember second field is the inserted one
    % So this is the anti region
    figure(1)
    mleplos(thhat1.*scl1,th1,[],[],[],[],[],p,sprintf('anti %s',p.mask),[])
    figure(1)
    figdisp(sprintf('%s_2a',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    clf
    figure(2)
    figdisp(sprintf('%s_2b',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    clf
    % An this is the selected region
    figure(1)
    mleplos(thhat2.*scl2,th2,[],[],[],[],[],p,p.mask,[])
    figure(1)
    figdisp(sprintf('%s_3a',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)
    figure(2)
    figdisp(sprintf('%s_3b',pref(sprintf('%s_%s.mat','MASKIT',fname))),[],[],2)

    % With PARFOR none of the once-used are available out of the loop but in
    % this demo2 you don't want any output anyway, so put in empties
    [v,cr,I,w,vw]=deal(NaN);
end

% Variable output
varns={v,cr,I,w,vw};
varargout=varns(1:nargout);

function plotit(v,p,cr,th)
% Should use th(1) instead of halverange, shouldn't I...
imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],halverange(v,80,NaN)); axis ij image 
hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
if ~isempty(th)
    t=title(sprintf('%i x %i | %s = [%g %g %gx]',p.NyNx(1),p.NyNx(2),'\theta',...
                    th./[1 1 sqrt(prod(p.dydx))]));
    movev(t,-p.NyNx(1)/20)
end
