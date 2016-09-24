function s=v2s(v)
% s=v2s(v)
%
% Vector-to-square matrix transform
%
% Last modified by fjsimons-at-alum.mit.edu, 06/17/2010

s=reshape(v,sqrt(length(v(:))),[]);
