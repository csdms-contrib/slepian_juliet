function varargout=mlexplos(thhats,covth,flabs)
% ah=MLEXPLOS(thhats,covth,flabs)
%
% Makes cross-plots of variables with their error ellipses, centered on the mean
% of a set or on individual values. If there are enough covariances given for
% every one of M sets of N variables, plots the entire population with their
% individual error ellipses. If one covariance is given, plots that for the
% population. In that case and in addition it calculates and plots one from the
% one set of M experiments of N variables.
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
% % Bivariate collective example
% mlexplos; axis square; axis([-3 3 -3 3]); grid on
% mlexplos(randx([1 1 0.5],128),[1 1 -0.5]); axis square; axis([-3 3 -3 3]); grid on
% % Trivariate collective example
% cv=[1 1 1 0.2 -0.6 0.6]; M=100; ah=mlexplos(randx(cv,M),cv,{'X','Y','Z'}); shrink(ah,1,2)
% % Trivariate set example
% thhats=[1 2 3 ; 4 5 6 ; 7 8 9 ; 10 11 12];
% cv=[1 1 1 0.2 -0.6 0.6 ; 1 1 1 -0.2 0.6 0.3 ; 1 1 1 -0.2 -0.6 0.9 ; 1 1 1 0.8 -0.7 0.6];
% mlexplos(thhats,cv,{'\sigma^2','\rho','\nu'})
%
% SEE ALSO:
%
% RANDX
%
% Last modified by fjsimons-at-alum.mit.edu, 07/20/2026

defval('covth',[1 1 0.5])
defval('thhats',randx(covth(1,:),100))
defval('flabs',{'X' 'Y'})

% The number of parameters
np=size(thhats,2);

% Number of pairwise combination plots
if np>2
    pcomb=nchoosek(1:np,2);
else
    pcomb=1;
end

% Try to make some arrangement of windows
clf
try
    [ah,ha]=krijetem(subnum(size(pcomb,1)/3,3));
catch
    try
        [ah,ha]=krijetem(subnum(size(pcomb,1)/2,3));
    catch
        ah=gca; ha=ah;
    end
end

% Single set with common covariance or collective with individual covariances
ss=prod(size(covth))==size(thhats,2)*(size(thhats,2)+1)/2;

% The below ripped off from MLEPLOS, which later should use MLEXPLOS
% Number of times the standard deviation for scale truncation
pstats=[-1 1];

for ind=1:np
  % The empirical means and standard deviations of the estimates
  if ss==1
      % Need to calculate the means and standard deviation of the collective
      mobss(ind)=nanmean(thhats(:,ind));
      sobss(ind)=nanstd(thhats(:,ind));
  else
      % You were given the center point as the single estimate each time
      mobss=thhats;
      for jnd=1:size(thhats,1)
          % Use the supplied standard deviation for each of the estimates
          sobss(jnd,ind)=sqrt(covth(jnd,ind));
      end
  end
end

% SCALING? 

% Plot error ellipses without using ERROR_ELLIPSE
% https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
defval('cl',0.68)

for ind=1:size(pcomb,1)
    axes(ah(ind))
    % Find the pairwise combination for the cross-plots
    p1=pcomb(ind,1); try p2=pcomb(ind,2); catch p2=2; end
    % Find the covariance matrix slice for this variable pair
    znp=zeros(1,np); znp([p1 p2])=1;
    
    % APPLY SCALING?

    % PLOT CROSSES WITH MEANS AND STANDARD DEVIATIONS?

    % The parameter estimates
    p(ind)=plot(thhats(:,p1),thhats(:,p2),'o');

    hold on
    % OBSERVED MEANS/ESTIMATES AND OBSERVED/CALCULATED STANDARD DEVIATIONS
    if ss==1
        m(ind)=plot(mobss(p1),mobss(p2),'v');
    end
    if ss==1; jmax=1; else jmax=size(thhats,1); end
    for jnd=1:jmax
        % Horizontal crosshair on the collective
        o1(jnd,ind)=plot(mobss(jnd,p1)+pstats*sobss(jnd,p1),[mobss(jnd,p2) mobss(jnd,p2)]);
        % Vertical crosshair on the collective
        o2(jnd,ind)=plot([mobss(jnd,p1) mobss(jnd,p1)],mobss(jnd,p2)+pstats*sobss(jnd,p2));
    end

    % FANCY TICKS AND LABELING

    % Estimates
    set(p(ind),'MarkerFaceColor',grey,'MarkerEdgeColor',grey,'MarkerSize',2)
    if ss==1
        % Means
        set(m(ind),'MarkerFaceColor',grey,'MarkerEdgeColor',grey,'MarkerSize',4)
    end
    % Crosshairs
    set([o1(ind) o2(ind)],'LineWidth',1,'Color',grey)
    % Ordering
    if ss==1
        top(m(ind),ah(ind))
    end
    for jn=1:size(thhats,1)
        bottom(o1(jnd,ind),ah(ind))
        bottom(o2(jnd,ind),ah(ind))
    end

    % Labels
    xl2(ind)=xlabel(flabs{p1});
    ylabel(flabs{p2})

    if ss==1
        % OBSERVED pairwise error ellipse for the collective
        % https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
        hold on
        disp('OBSERVED')
        ep(ind)=covell(cl,cov(thhats(:,[p1 p2])),thhats(:,[p1 p2]));
        
        % SUPPLIED pairwise error ellipse for the collective
        disp('CALCULATED')
        ec(ind)=covell(cl,matslice(trilosi(covth),znp),thhats(:,[p1 p2]));

        % Observed covariance ellipse
        set(ep(ind),'LineWidth',1.5,'Color',grey)
        % Supplied covariance ellips
        set(ec(ind),'LineWidth',0.5,'Color','k')
        
        % Send these ellipses to the back so the dots show on top
        bottom(ec(ind),ah(ind))
        bottom(ep(ind),ah(ind))
    else
        % They are all different estimates with their uncertainty calculations
        % SUPPLIED pairwise error ellipse for the individuals
        disp('CALCULATED')
        if ss==1
            ec(ind)=covell(cl,...
                           matslice(trilosi(covth(jnd,:)),znp),thhats(:,[p1 p2]));
        else
            for jnd=1:size(thhats,1)
                ec(jnd,ind)=covell(cl,...
                                   matslice(trilosi(covth(jnd,:)),znp),thhats(jnd,[p1 p2]));
            end
        end
    end

    if ss==1
        % These things often normalized so this would be appropriate
        axis image; axis([-3 3 -3 3]); grid on
    end
        
    % Send the grid lines all the way to the back for FANCY TICKS AND LABELING
    hold off
end

% Cosmetix
longticks(ah)

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
if size(thhats,1)>1
    % Calculate the one center
    mobs=nanmean(thhats,1);
else
    % The inputs are the individual centers
    mobs=thhats;
end

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



