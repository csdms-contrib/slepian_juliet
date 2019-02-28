function [pf,pv,pr,alfa]=normtest(stat,mn,vr,alfa)
% [pf,pv,pr,alfa]=NORMTEST(stat,mn,vr,alfa)
%
% Test whether a particular value is derived from a hypothesized NORMAL
% (GAUSSIAN) distribution with known (postulated) mean and variance, reject
% if the probability of more extreme than observed values is very unlikely,
% e.g. below (1-alfa)x100 per cent. This is a TWO-SIDED test, where the
% algorithm depends on the symmetry of the distribution. It is a ONE-SAMPLE
% test for every ELEMENT of the array offered (i.e. quite literally ONE
% sample, of one point, each), unlike Matlab's own ZTEST.
%
% INPUT: 
%
% stat      The values being tested individually, could be a vector or a matrix
% mn        The expected value of the normal distribution under the null [0]
% vr        The variance of the  normal distribution under the null [1]
% alfa      The significance level for a (1-alfa)x100 confidence [0.05]
%
% OUTPUT:
%
% pf        0 if it PASSES, i.e. the statistic is derived from the null
%           1 if it FAILS, i.e. the statistic is REJECTED to be from null 
% pv        The 'p-value' of having an even more extreme value of the test
% pr        The percentage of the test values that are being rejected,
%           useful if you give a whole vector of stat values
% alfa      The significance level [0.05] for a (1-alfa)x100 confidence
% 
% EXAMPLE:
%
% mu=randn; vr=rand; al=randi(100)/100; mn=randi(1e3); [a,b,c]=normtest(mu+sqrt(vr)*randn(mn,1),mu,vr,al);
%
% SEE ALSO: 
%
% ZTEST, VARTEST
%
% Last modified by fjsimons-at-alum.mit.edu, 07/03/2018

% Default confidence limit
defval('alfa',0.05)

% ONE-SAMPLE (single value!) TEST ON THE MEAN (applied on every element)
% This is the probability of a single test value even more extremely
% removed, in the absolute sense, from the expectation, given a known
% standard deviation 
pv=1-2*abs(normcdf(stat,mn,sqrt(vr))-1/2);

% This is the test, reject gives a 1, as in, could arise by chance
pf=pv<alfa;

% This is the report on the rejections
pr=sum(pf)/length(pf)*100;

% Provide an uplifting message
if length(stat)>1
  disp(sprintf('\n%s %i : %i%% at the %i%% confidence level rejected for m %8.3f v %8.3f',...
	       upper(mfilename),length(stat),round(pr),round(alfa*100),mn,vr)) 
end
