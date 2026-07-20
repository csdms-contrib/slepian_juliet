function varargout=mlexplos(thhats,covth,flabs)
% ah=MLEXPLOS(thhats,covth,flabs)
%
% Makes cross-plots of variables with their error ellipses, centered on the
% mean. If there are enough covariances given for every one of M sets of N
% variables, plots the entire population with their individual error ellipses.
% If no covariance is given, calculates one from the one set of M experiments of N variables.
%
% INPUT:
%
% thhats    A MxN dimensional set of N variables
% covth     A MpxNp rolled-up covariance matrix (matrix), calculated if unavailable
%               NP=N*(N+1)/2 unique entries in an NxN symmetric matrix, see TRILOS(I)
%               Mp is 1 if it pertains to the entire set of thhats
%               Mp is M if it pertains to the entries in the set of thhats individually
% flabs     Cell with label strings for each of the variables
%
% OUTPUT:
%
% ah        Axis handles
%
% EXAMPLE:
%
% mlexplos; axis square; axis([-3 3 -3 3]); grid on
% mlexplos(randx([1 1 1 0.2 -0.6 0.6]),[1 1 1 0.2 -0.6 0.6],{'X','Y','Z'})
%
% SEE ALSO:
%
% RANDX
%
% Last modified by fjsimons-at-alum.mit.edu, 07/20/2026

defval('covth',[1 1 0.5])
defval('thhats',randx(covth,100))
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
clf
keyboard
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
    % Horizontal crosshair
    o1(ind)=plot(mobss(p1)+pstats*sobss(p1),[mobss(p2) mobss(p2)]);
    % Vertical crosshair
    o2(ind)=plot([mobss(p1) mobss(p1)],mobss(p2)+pstats*sobss(p2));

    % FANCY TICKS AND LABELING

    % Estimates
    set(p(ind),'MarkerFaceColor',grey,'MarkerEdgeColor',grey,'MarkerSize',3)
    % Means
    set(m(ind),'MarkerFaceColor',grey,'MarkerEdgeColor',grey,'MarkerSize',4)
    % Crosshairs
    set([o1(ind) o2(ind)],'LineWidth',1,'Color',grey)
    % Ordering
    top(m(ind),ah(ind))
    bottom(o1(ind),ah(ind))
    bottom(o2(ind),ah(ind))

    % Labels
    xl2(ind)=xlabel(flabs{p1});
    ylabel(flabs{p2})

    % OBSERVED pairwise error ellipse
    % https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
    hold on
    disp('OBSERVED')
    ep(ind)=covell(cl,cov(thhats(:,[p1 p2])),thhats(:,[p1 p2]));
    
    % CALCULATED pairwise error ellipse
    znp=zeros(1,np); znp([p1 p2])=1;
    disp('CALCULATED')
    ec(ind)=covell(cl,matslice(covth,znp),thhats(:,[p1 p2]));

    % Observed covariance ellipse
    set(ep(ind),'LineWidth',1.5,'Color',grey)
    % Supplied covariance ellips
    set(ec(ind),'LineWidth',0.5,'Color','k')
    
    % Send these ellipses to the back so the dots show on top
    bottom(ec(ind),ah(ind))
    bottom(ep(ind),ah(ind))

    % Send the grid lines all the way to the back for FANCY TICKS AND LABELING

    hold off
end

% Optional outputs
varns={ah};
varargout=varns(1:nargout);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=covell(cl,covx,thhats)
% Plots bivariate covariance ellipse and reports on inside counts

% Check this for three variables
% s=-2*log(1-cl);
s=chi2inv(cl,2);
% Make circular coordinates
t=linspace(0,2*pi);
% Calculate center
mobs=nanmean(thhats,1);

% Compute the eigenvectors and eigenvalues of the covariance
[V,D]=eig(covx);
% Compute coordinates of the ellipse
a=sqrt(s)*V*sqrt(D)*[cos(t); sin(t)];
% Actually do the plotting
p=plot(a(1,:)+mobs(1),a(2,:)+mobs(2));

% Count the number of estimates outside the OBSERVED confidence interval
cli=round(sum(inpolygon(thhats(:,1),thhats(:,2),...
                        a(1,:)+mobs(1),a(2,:)+mobs(2)))/size(thhats,1)*100);
disp(sprintf('CL asked %g ; received %g ',cl*100,cli))



