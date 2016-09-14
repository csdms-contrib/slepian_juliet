function hermcheck(A)
% HERMCHECK(A)
%
% Check the Hermitian symmetry of a spectral matrix
%
% INPUT:
%
% A     The matrix
%
% SEE ALSO:
%
% KNUM2
%
% Last modified by fjsimons-at-alum.mit.edu, 05/04/2016

% Disregarding the first row and column if the size is even
Hh=A(1+~mod(size(A,1),2):end,1+~mod(size(A,2),2):end);
diferm(conj(fliplr(flipud(Hh)))-Hh,[],12-round(log10(norm(Hh))))
