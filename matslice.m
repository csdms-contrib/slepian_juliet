function B=matslice(A,keepdim,r)
% B=matslice(A,KEEPDIM,R)
%
% Remove row(s) and column(s) etc from an array
%
% INPUTS:
%
% A         the full matrix
% keepdim   the logical flags for kept dimensions parameters;
%           removal will happen for 0's
% r     -1  reverse the slice; introduce row(s) and col(s) of NaNs
%
% OUTPUTS:
%
% B     the sliced matrix
%
% SEE ALSO:
% 
% VARIANCEOFESTIMATES, MLEOSL
%
% EXAMPLE:
%
% A=rand(3); keepdim=[1 0 1];
% B=matslice(A,keepdim);
% disp(sprintf('Size check: %s',string(all(size(B)==[2,2]))))
% 
% A=rand(3,1);
% B=matslice(A,keepdim);
% disp(sprintf('Size check: %s',string(all(size(B)==[2,1]))))
%
% A=rand(1,3);
% B=matslice(A,keepdim);
% disp(sprintf('Size check: %s',string(all(size(B)==[1,2]))))
% C=B;
% B=matslice(C,keepdim,-1);
% disp(sprintf('Size check: %s',string(all(size(B)==size(A)))))
%
% matslice(matslice(randn(3,3),[1 0 1]),[1 0 1],-1)
%
% Last modified by fjsimons-at-princeton.edu, 04/08/2025
% Last modified by olwalbert-at-princeton.edu, 04/08/2025

% Create a cell array to match the dimensions of A with logicals as elements 
% for indexing
defval('r',1)

idx=repmat({logical(keepdim)},sum(size(A)>1),1);

if r==1
  % How to 'index' rows and columns etc all at the same time 
  B=A(idx{:});
elseif r==-1
  np=size(keepdim,2);
  dm=sum(size(A)>1);
  B=squeeze(reshape(zeros(1,np.^dm),[1 repelem(np,dm)]));
  % Would be nice to get NaNs for blanks, but TRILOS doesn't like it
  % B(:)=deal(NaN);
  B(idx{:})=A;
end
