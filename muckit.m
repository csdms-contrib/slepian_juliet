function varargout=muckit(v,p,scl)
% [v,cr,I]=MUCKIT(v,p,scl)
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
%
% OUTPUT:
%
% v       The masked output matrix, unrolled into a vector
% cr      The colum,row index coordinates of the bounding curve
% I       The mask, unrolled into a vector, note, this is an "anti-mask"
% scl     A scale [0 to 1] expressing the non-missing data fraction%
%
% EXAMPLE:
%
% [~,~,I]=muckit('demo1','random',0.8); % To get a quick muck to look at
% muckit('demo2'); % To do a serious of mock muck inversions
%
% Last modified by fjsimons-at-alum.mit.edu, 10/13/2023

% The default is the demo, for once
defval('v','demo1')

if ~isstr(v)
    % Get the muck by name? Let's play "mikado"?
    if isstr(p.mask)
        % Scale indicates survivorship (in symmetry with MASKIT)
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
    % Apply the mask to the first field
    v(~I(:))=0;
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
    [Hm,cr,I]=muckit(Hx,p,scl);
    clf
    subplot(121); plotit(Hx,p,cr,th)
    subplot(122); plotit(Hm,p,cr,[])
elseif strcmp(v,'demo2')
    % Capture the second input
    defval('p','random')
    defp=p; clear p
    % Now proceed with a fresh copy
    p.mask=defp;
    % Simulate using invariant embedding, no taper
    p.quart=0; p.blurs=Inf; p.kiso=NaN; clc;
    % Something manageable without overdoing it
    p.NyNx=[188 233]+randi(20,[1 2]);
    % Something larger without overdoing it, check weirdness
    p.NyNx=[512 512]/8;

    N=3;
    for index=1:N
        clc; disp(sprintf('\n Simulating the field \n'))

        % Simulate the field with defaults from SIMULOSL
	% and collect the expected periodogram...
        [Hx,th1,p,~,~,Sb]=simulosl([],p,1);

        % How much should the surviving area occupy?
        scl=[];
        % Now do the masking and the merging
        [Hm,cr,I,scl]=muckit(Hx,p,scl);

        % Make a visual for good measure
        clf
        ah(1)=subplot(121); plotit(Hx,p,[],th1)
	xlabel('Original field')
        ah(2)=subplot(122); plotit(Hm,p,[],th1); 
	xlabel(sprintf('%i%% speckled field',round((1-scl)*100)))

        % So HG is the mixed field
        v=Hm;

        % As appropriate you'll force the use of BLUROSY in MATERNOSP in LOGLIOS
        p.blurs=-1;
        % Do all the tests or not
        xver=0;

        % Recover the parameters of the full original field without any masking
        p.taper=0;
        % Make a close initial guess?
        thini=th1+(-1).^randi(2,[1 3]).*th1/1000;
        pause(5); clc; disp(sprintf('\n Estimating first whole field \n'))
        [thhat1(index,:),~,~,scl1(index,:)]=mleosl(Hx,thini,p,[],[],[],xver);
	% Explicitly check the likelihood?
	disp(sprintf('\nLast check on likelihood at the initial guess\n'))
	logliosl(knums(p),thini./scl1(index,:),scl1(index,:),p,tospec(Hx(:),p)/(2*pi));
	disp(sprintf('\nLast check on likelihood at the truth\n'))
	logliosl(knums(p),th1./scl1(index,:),scl1(index,:),p,tospec(Hx(:),p)/(2*pi));
	% --> inside protect by looking at momx?? That's where it goes wrong

        % Now recover the parameters of the speckled field
        p.taper=I;
        % Make a close initial guess?
        thini=th1+(-1).^randi(2,[1 3]).*th1/1000;
        pause(5); clc; disp(sprintf('\n Estimating first speckled field \n'))
        [thhat2(index,:),~,~,scl2(index,:)]=mleosl(Hm,thini,p,[],[],[],xver);
	% Explicitly check the likelihood?
	disp(sprintf('\nLast check on likelihood at the initial guess\n'))
	logliosl(knums(p),thini./scl1(index,:),scl1(index,:),p,...
		 tospec(p.taper(:).*Hx(:),p)/(2*pi)/sqrt(sum(p.taper(:).^2))*sqrt(prod(p.NyNx)),1);
	disp(sprintf('\nLast check on likelihood at the truth\n'))
	logliosl(knums(p),th1./scl1(index,:),scl1(index,:),p,...
		 tospec(p.taper(:).*Hx(:),p)/(2*pi)/sqrt(sum(p.taper(:).^2))*sqrt(prod(p.NyNx)),1);
    end
    % And now look at the statistics of the recovery
    disp(sprintf('\nWhole | Speckled\n'))
    disp(sprintf('%8.0f %5.2f %6.0f  %8.0f %5.2f %6.0f\n',[thhat1.*scl1 thhat2.*scl2]'))
    
    % Then plot these things using MLEPLOS
    mleplos(thhat1,th1,[],[],[],[],[],p,'MUCKIT-demo2')
end

% Variable output
varns={v,cr,I,scl};
varargout=varns(1:nargout);

function plotit(v,p,cr,th)
% Should use th(1) instead of halverange, shouldn't I...
imagefnan([1 1],p.NyNx([2 1]),v2s(v,p),[],halverange(v,80,NaN)); axis ij image 
%hold on; twoplot(cr,'Color','k'); hold off; longticks; grid on
xticks([1 round(p.NyNx(2)/2) p.NyNx(2)]); yticks([1 round(p.NyNx(1)/2) p.NyNx(1)])
if ~isempty(th)
    t=title(sprintf('%i x %i | %s = [%g %g %gx]',p.NyNx(1),p.NyNx(2),'\theta',...
                    th./[1 1 sqrt(prod(p.dydx))]));
    movev(t,-p.NyNx(1)/20)
end
