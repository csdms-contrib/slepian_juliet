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
% maskit('demo3') % Plotting results of previous demo
%
% Last modified by fjsimons-at-alum.mit.edu, 05/31/2024

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

    if xver==1
        pz=2;
    else
        pz=0.1;
    end

    N=10;
    for index=1:N
        clc; disp(sprintf('\n Simulating all fields \n'))

        % Define the parameters for the two fields
        th1=[1.5 1.25 30000];
        th2=[1.9 1.75 10000];

        % Simulate first field
        [Hx,th1,p,k,Hk,Sb1]=simulosl(th1,p,1);
        % Simulate second field
        [Gx,th2,p,k,Hk,Sb2]=simulosl(th2,p,1);

        if xver==1
            % Very explicit comparison like BLUROSY
            p.blurs=-1; Sbar1=blurosy(th1,p,1); 
            p.blurs=-1; Sbar2=blurosy(th2,p,1);
            % Are any of them bad?
            if max(Sb1./Sbar1)>20;keyboard ;end
            if max(Sb2./Sbar2)>20;keyboard ;end
            % Reset to the original value
            p.blurs=Inf;
        end
        
        % How much should the masked area occupy?
        scl=[];
        % Now do the masking and the merging
        [Hm,cr,I,Gm,HG]=maskit(Hx,p,scl,Gx);

        % Save data for subsequent visualization 
        % save ICFrance p HG I cr th1 th2

        % Make a visual for good measure
        clf
        ah(1)=subplot(221); plotit(Hx,p,cr,th1)
        ah(3)=subplot(223); plotit(Gx,p,cr,th2)
        ah(4)=subplot(224); plotit(HG,p,cr,[])
        movev(ah(4),0.25)
        % So HG is the mixed field
        v=Hm; w=Gm; vw=HG;

        % Without the boundary the eye is much less clear!
        % imagesc(v2s(HG,p)); axis image

        % As appropriate you'll force the use of BLUROSY in MATERNOSP in LOGLIOSL in MLEOSL
        p.blurs=-1;

        % Make a close initial guess?
        thini1=th1+(-1).^randi(2,[1 3]).*th1/1000;
        thini2=th2+(-1).^randi(2,[1 3]).*th2/1000;

        % Recover the parameters of the full original fields without any masking
        p.taper=0;
        pause(pz); clc; disp(sprintf('\n Estimating first whole field \n'))
        [thhat4(index,:),~,~,scl4(index,:)]=mleosl(Hx,thini1,p,[],[],[],xver);
        pause(pz); clc ; disp(sprintf('\n Estimating second whole field \n'))
        [thhat5(index,:),~,~,scl5(index,:)]=mleosl(Gx,thini2,p,[],[],[],xver);

        % Now recover the parameters of the mixed field but only in the region I or ~I
        % Perform the optimization on the insert which should look like the second field
        p.taper=I;
        pause(pz); clc; disp(sprintf('\n Estimating first partial field \n'))
        [thhat1(index,:),~,~,scl1(index,:)]=mleosl(HG,thini2,p,[],[],[],xver);
        
        % Take a look inside LOGLIOS that the order of magnitude is good.
        
        % Now recover the parameters of the complement which should look like the first field
        p.taper=~I;
        pause(pz); clc; disp(sprintf('\n Estimating second partial field \n'))
        [thhat2(index,:),~,~,scl2(index,:)]=mleosl(HG,thini1,p,[],[],[],xver);

        % Now recover the parameters of the merged field without knowing of the partition
        % p.taper=0;
        % pause(pz); clc; disp(sprintf('\n Estimating whole mixed field \n'))
        % [thhat3(index,:),~,~,scl3(index,:)]=mleosl(HG,[],p,[],[],[],xver);

        % Pause so you can watch live (big pause if xver==0)
        pause(pz)

        % Reset to the original values
        p.blurs=Inf; p.taper=0;
    end
    % And now look at the statistics of the recovery
    disp(sprintf('\nFirst partial | Second partial\n'))
    disp(sprintf('%8.2f %5.2f %6.0f  %8.2f %5.2f %6.0f\n',[thhat1.*scl1 thhat2.*scl2]'))
    disp(sprintf('\nFirst whole | Second whole\n'))
    disp(sprintf('%8.2f %5.2f %6.0f  %8.2f %5.2f %6.0f\n',[thhat4.*scl4 thhat5.*scl5]'))
    %disp(sprintf('\nMixed whole\n'))
    %disp(sprintf('%8.0f %5.2f %6.0f\n',[thhat3.*scl3]'))
    
    % Then plot these things using MLEPLOS
    mleplos(thhat1.*scl1,th1,[],[],[],[],[],p,[],[])
    mleplos(thhat2.*scl2,th2,[],[],[],[],[],p,[],[])

    % Save some output for later?
    keyboard
elseif strcmp(v,'demo3')
    % Saved by hand and plot at the end
    maskitdemo2_10012023
    th1=[1e+06  2.5  20000];
    th2=[1e+06  2.5  30000];
    p.NyNx=[205 252];
    p.dydx=[10000 10000];
    p.blurs=Inf; p.kiso=NaN;
    % The first partial contains the second field
    mleplos(trimit(FirstPartialSecondPartial(:,1:3),97),th2,[],[],[],[],[],p,'MASKIT-10012023-I')
    figure(1)
    figdisp('maskitdemo2_10012023_Ia',[],[],2)    
    figure(2)
    figdisp('maskitdemo2_10012023_Ib',[],[],2)
    % The second partial contains the first field
    mleplos(trimit(FirstPartialSecondPartial(:,4:6),97),th1,[],[],[],[],[],p,'MASKIT-10012023-II')
    figure(1)
    figdisp('maskitdemo2_10012023_IIa',[],[],2)   
    figure(2)
    figdisp('maskitdemo2_10012023_IIb',[],[],2)
    
    % The first whole
    mleplos(trimit(FirstWholeSecondWhole(:,1:3),97),th1,[],[],[],[],[],p,'MASKIT-10012023-0')
    figure(1)
    figdisp('maskitdemo2_10012023_0a',[],[],2)    
    figure(2)
    figdisp('maskitdemo2_10012023_0b',[],[],2)
    % The second partial contains the first field
    mleplos(trimit(FirstWholeSecondWhole(:,4:6),97),th2,[],[],[],[],[],p,'MASKIT-10012023-00')
    figure(1)
    figdisp('maskitdemo2_10012023_00a',[],[],2)   
    figure(2)
    figdisp('maskitdemo2_10012023_00b',[],[],2)
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
