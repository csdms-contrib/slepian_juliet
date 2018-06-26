function s=v2s(v,varargin)
% s=v2s(v,varargin)
%
% Vector-to-square/rectangle matrix transform
%
% INPUT:
%
% v       A vector that is the unwrapping of a matrix
% p       Optional parameter structure as from the MLEOSL suite, at least:
%         NyNx  number of samples in the y and x directions
%
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2018

if nargin==1
  s=reshape(v(:),sqrt(length(v(:))),[]);
else
  s=reshape(v(:),varargin{1}.NyNx);
end
