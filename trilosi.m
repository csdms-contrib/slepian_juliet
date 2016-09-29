function m=trilosi(v)
% m=trilosi(v)
% 
% Returns a full symmetric matrix from a version slimmed down by TRILOS
%
% INPUT:
%
% v      The elements that matter, linearly unwrapped
% np     The number of diagonal elements
%
% OUTPUT:
%
% m      The full matrix
% 
% SEE ALSO:
%
% TRILOS
%
% Last modified by fjsimons-at-alum.mit.edu, 10/25/2014

% How many diagonal elements?
np=(-1+sqrt(1+4*2*length(v(:))))/2;

% Empty bucket
m=zeros(np,np);

% Fill the lower half of them up
m(nonzeros(triu(reshape(1:np^2,np,np)')'))=v;

% Symmetrize
m=[tril(m)'+tril(m)-diag(diag(m))];
