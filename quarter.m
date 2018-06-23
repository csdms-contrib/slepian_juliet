function M=quarter(M)
% M=QUARTER(M)
%
% Extracts the central 'quarter' of a 'divisibly' sized matrix
%
% INPUT:
%
% M   The original matrix
%
% OUTPUT:
%
% M   The new matrix
%
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2018

[m,n]=size(M);

warning off MATLAB:colon:nonIntegerIndex 
M=M(m/4+1:3*m/4,n/4+1:3*n/4);
warning on MATLAB:colon:nonIntegerIndex 
