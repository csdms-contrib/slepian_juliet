function mlexplos(thhats,covth)
% MLEXPLOS(thhats,covth)
%
% Makes cross-plots of variables with their error ellipses. If there are enough
% covariances given for every one of M sets of N variables, plots the entire
% population with their individual error ellipses. If no covariance is given,
% calculates one from the one set of M experiments of N variables.
%
% INPUT:
%
% thhats    A MxN dimensional set of N variables
% covth     A MpxNp rolled-up covariance matrix (matrix), calculated if unavailable
%               NP=N*(N+1)/2 unique entries in an NxN symmetric matrix, see TRILOS(I)
%               Mp is 1 if it pertains to the entire set of thhats
%               Mp is M if it pertains to the entries in the set of thhats individually
%
% EXAMPLE:
%
% mlexplos; axis square; axis([-3 3 -3 3]); grid on
%
% SEE ALSO:
%
% RANDX
%
% Last modified by fjsimons-at-alum.mit.edu, 07/20/2026

defval('thhats',randx)
defval('covth',[1 1 0.5])
defval('flabs',{'X' 'Y'})

% The number of parameters
np=size(thhats,2);

% Number of pairwise combination plots
if np>2
    pcomb=nchoosek(1:np,2);
else
    pcomb=1;
end

% Rearrange into proper symmetric full form
covth=trilosi(covth);

% Try to make some arrangement of windows
try
    [ah,ha]=krijetem(subnum(size(pcomb,1)/3,3));
catch
    try
        [ah,ha]=krijetem(subnum(size(pcomb,1)/2,3));
    catch
        ah=gca; ha=ah;
    end
end

% The below ripped off from MLEPLOS, which later should use MLEXPLOS
% Number of times the standard deviation for scale truncation
nstats=[-3:3]; fax=3;
pstats=[-1 1];

for ind=1:np
  % The empirical means and standard deviations of the estimates
  % Collect them all
  mobss(ind)=nanmean(thhats(:,ind));
  sobss(ind)=nanstd(thhats(:,ind));
end

% SCALING? 

% Plot error ellipses without using ERROR_ELLIPSE
% https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
defval('cl',0.68)
% Check this for three variables
% s=-2*log(1-cl);
s=chi2inv(cl,2);

t=linspace(0,2*pi);

for ind=1:size(pcomb,1)
    axes(ah(ind))
    % Find the pairwise combinations for the cross-plot convention:

    p1=pcomb(ind,1); try p2=pcomb(ind,2); catch p2=2; end

    % APPLY SCALING?

    % PLOT CROSSES WITH MEANS AND STANDARD DEVIATIONS?

    % The parameter estimates
    p(ind)=plot(thhats(:,p1),thhats(:,p2),'o');

    hold on
    % OBSERVED MEANS AND OBSERVED STANDARD DEVIATIONS
    m(ind)=plot(mobss(p1),mobss(p2),'v');
    o1(ind)=plot(mobss(p1)+pstats*sobss(p1),...
		 [mobss(p2) mobss(p2)],'LineWidth',1);
    o2(ind)=plot([mobss(p1) mobss(p1)],...
		 mobss(p2)+pstats*sobss(p2),'LineWidth',1);

    % FANCY TICKS AND LABELING

    % Color mix
    cmix=[0 0 0]; 
    set([p(ind) m(ind)],'MarkerFaceColor',cmix,'MarkerEdgeColor',cmix,'MarkerSize',2)
    % Dull colors
    set([p(ind) m(ind)],'MarkerFaceColor',grey,'MarkerEdgeColor',grey,'MarkerSize',2)

    % Cosmetix
    % Emphasize the OBSERVED means and standard deviations
    set([o1(ind) o2(ind)],'LineWidth',1,'Color',grey(3.5))
    set( m(ind)           ,'MarkerFaceColor',grey(3.5),'MarkerEdgeColor',grey(3.5),'MarkerSize',4)
    top(m(ind),ah(ind))
    bottom(o1(ind),ah(ind))
    bottom(o2(ind),ah(ind))

    % Labels
    xl2(ind)=xlabel(flabs{p1});
    ylabel(flabs{p2})

    % OBSERVED
    % Plot pairwise error ellipses
    % https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
    hold on
    % Compute the eigenvectors and eigenvalues of the covariance
    % Think of the Schur complement? How does that relate?
    [V,D]=eig(cov(thhats(:,[p1 p2])));
    a=sqrt(s)*V*sqrt(D)*[cos(t); sin(t)];
    ep(ind)=plot(a(1,:)+mobss(p1),a(2,:)+mobss(p2));

    % Count the number of estimates outside the confidence interval
    cli=round(sum(inpolygon(thhats(:,p1),thhats(:,p2),...
                            a(1,:)+mobss(p1),a(2,:)+mobss(p2)))/size(thhats,1)*100);
    disp(sprintf('CL asked %g ; received %g ',cl*100,cli))

    
    % And for the calculated covariance, too
    znp=zeros(1,np);
    znp([p1 p2])=1;
    [V,D]=eig(matslice(covth,znp));
    % Get rid of the a11 line in MLEPLOS
    a=sqrt(s)*V*sqrt(D)*[cos(t); sin(t)];
    ec(ind)=plot(a(1,:)+mobss(p1),a(2,:)+mobss(p2));

    % Dull colors
    set(ep(ind),'LineWidth',1.5,'Color',grey(3.5))
    set(ec(ind),'LineWidth',1,'Color','k')
    % Send these ellipses to the back so the dots show on top
    bottom(ec(ind),ah(ind))
    bottom(ep(ind),ah(ind))

    % Send the grid lines all the way to the back for FANCY TICKS AND LABELING

    keyboard

    hold off
end





