function v=trilos(m)
% v=trilos(m)
% 
% Slims down a symmetric matrix for outputting covariances, etc.
%
% INPUT:
%
% m       A full-form symmetric matrix (and checked to be so)
%
% OUTPUT:
%
% v      The elements that matter, linearly unwrapped
%
% SEE ALSO:
%
% TRILOSI
%
% Last modified by fjsimons-at-alum.mit.edu, 02/09/2015

% Check symmetry
difer(m-m',[],[],NaN)

% Keep the lower-triangular portion
m=tril(m); 

% Linearly unwrap it
v=m(~~m);
