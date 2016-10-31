function v=trilos(m)
% v=trilos(m)
% 
% Slims down a symmetric (e.g. covariance, Fisher or Hessian) matrix by
% first going down the diagonal and then across the remaining columns
% down the rows.
%
% INPUT:
%
% m       A full-form symmetric matrix
%
% OUTPUT:
%
% v      The elements that matter, in a column vector
%
% SEE ALSO:
%
% TRILOSI, NCHOOSEK
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2016

% Do'nt check symmetry as we will "fake-use" TRILOS for TRILOSI
% difer(m-m',[],[],NaN)

% Reorganize
dm=diag(m);
m=tril(m)-diag(dm); 

% Linearly unwrap it
v=[dm ; m(~~m)];
