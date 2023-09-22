function varargout=maskit(v,p,scl,w)
% [v,cr,I,w,vw]=MASKIT(v,p,scl,w)
%
% INPUT:
%
% v       A vector that is the unwrapping of a matrix
% p       A single parameter structure with AT LEAST this:
%           NyNx  number of samples in the y and x directions
%           mask an index matrix with a mask, OR:
%                 a region name that will be scaled to fit
% scl     A scale between 0 and 1 with the occupying center
%         fraction of the domain enclosed by the input curve; due
%         to strict inequality at maximum always leaves one pixel
%         rim; and not being used if an actual mask is input
% w       A second vector that is the unwrapping of a second matrix
%
% OUTPUT:
%
% v       The masked output matrix, unrolled into a vector
% cr      The colum,row index coordinates of the masking curve
% I       The mask, unrolled into a vector, note, this is an "anti-taper"
% w       The second masked output matrix, unrolled into a vector
% vw      The merged field, second field w sucked into first field v
%
% EXAMPLE:
%
% [~,~,I]=maskit('demo1'); % To get a quick mask to look at
% maskit('demo1') % A geographical region masking a single field
% maskit('demo1','england') % 'amazon', 'orinoco', for geographical variability
% maskit('demo2') % A geographical region merging two fields
% maskit('demo2','england') % 'amazon', 'orinoco', for geographical variability
%
% Last modified by fjsimons-at-alum.mit.edu, 09/22/2023

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
        w(~I(:))=NaN;
        wv=NaN;
        if nargout>4
            vw=nan(size(v));
            vw(I(:)) =w(I(:));
            vw(~I(:))=v(~I(:));
        end
    end
    % This is the order, dudes
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
    % Simulate using invariant embedding, no taper
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc;
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);

    N=3;
    for index=1:N
        % Simulate first field with defaults from simulosl
        [Hx,th,p]=simulosl([],p,1);
        % Simulate second field by making a small change
        th2=th; th2(2)=2.5; th2(3)=40000;
        [Gx,th2,p]=simulosl(th2,p,1);

        % Now do the masking and the merging
        [Hm,cr,I,Gm,HG]=maskit(Hx,p,[],Gx);
        % No merging? [Hm,cr,I,Gm,HG]=maskit(Hx,p,[],Hx);
        clf
        ah(1)=subplot(221); plotit(Hx,p,cr,th)
        ah(3)=subplot(223); plotit(Gx,p,cr,th2)
        ah(4)=subplot(224); plotit(HG,p,cr,[])
        movev(ah(4),0.25)
        % So HG is the mixed field
        v=Hm; w=Gm; vw=HG;

        % Without the boundary the eye is much less clear!
        % imagesc(v2s(HG,p)); axis image

        warning('Not ready yet')
        try
        % As appropriate you'll force the use of BLUROSY in MATERNOSP in LOGLIOS
        p.blurs=-1;

        % Now recover the parameters of HG but only in the region I or ~I
        % Perform the optimization on the insert
        p.taper=I;
        [thhat(index,:),~,~,scl(index,:)]=mleosl(HG,[],p,[],[],[],[]);

        % Take a look inside LOGLIOS that the order of magnitude is good. Try for a close initial guess

        % Now recover the parameters of the complement
        p.taper=~I;
        [thhat2(index,:),~,~,scl2(index,:)]=mleosl(HG,[],p,[],[],[],[]);
        end
    end

    % And now look at the statistics of the recovery
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
    t=title(sprintf('%i x %i | %s = [%g %g %g]',p.NyNx(1),p.NyNx(2),'\theta',...
                    th./[1 1 sqrt(prod(p.dydx))]));
    movev(t,-p.NyNx(1)/20)
end
