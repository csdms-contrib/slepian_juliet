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
% mlexplos; axis square; axis([-3 3 -3 3])
%
% SEE ALSO:
%
% RANDX
%
% Last modified by fjsimons-at-alum.mit.edu, 07/20/2026

defval('thhats',randx)
defval('covth',[1 1 0.5])

% The number of parameters
np=size(thhats,2);

% Number of pairwise combination plots
if np>2
    pcomb=nchoosek(1:np,2);
else
    pcomb=1;
end

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
end





