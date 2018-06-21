function m=trilosi(v)
% m=trilosi(v)
% 
% Returns a full symmetric matrix from a version slimmed down by TRILOS
%
% INPUT:
%
% v      The elements that matter, in a column vector
%
% OUTPUT:
%
% m      The full matrix
% 
% SEE ALSO:
%
% TRILOS, NCHOOSEK
%
% Last modified by fjsimons-at-alum.mit.edu, 06/20/2018

% Figure out the full dimension of the matrix
np=(-1+sqrt(1+4*2*length(v(:))))/2;

% Empty bucket
m=zeros(np,np);

% Fill the lower half of them up
m(trilos(reshape(1:np^2,np,np)))=v;

% Symmetrize
m=[m+m'-diag(diag(m))];
