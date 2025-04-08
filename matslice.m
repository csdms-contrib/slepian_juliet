function B=ifinvslc(A,ifinv,r)
% B=IFINVSLC(A,IFINV,R)
%
% Slices the input matrix to remove the row(s) and column(s) that pertain to a
% fixed (not inverted) value. 
%
% INPUTS:
%
% A     the full matrix
% ifinv the logical flags for inverted parameters; slicing will occur for 0's
% r     -1 reverse the slice; introduce row(s) and col(s) of NaNs
%
% OUTPUTS:
%
% B     the sliced matrix
%
% SEE ALSO:
% 
% varianceofestimates, mleosl, matslice
%
% EXAMPLE:
%
% A=rand(3); ifinv=[1 0 1];
% B=ifinvslc(A,ifinv);
% disp(sprintf('Size check: %s',string(all(size(B)==[2,2]))))
% 
% A=rand(3,1);
% B=ifinvslc(A,ifinv);
% disp(sprintf('Size check: %s',string(all(size(B)==[2,1]))))
%
% A=rand(1,3);
% B=ifinvslc(A,ifinv);
% disp(sprintf('Size check: %s',string(all(size(B)==[1,2]))))
% C=B;
% B=ifinvslc(C,ifinv,-1);
% disp(sprintf('Size check: %s',string(all(size(B)==size(A)))))
%
% Last modified by olwalbert-at-princeton.edu, 01/28/2025

% Create a cell array to match the dimensions of A with logicals as elements 
% for indexing
defval('r',1)

idx=repmat({logical(ifinv)},sum(size(A)>1),1);

if r==1
  % How to 'index' rows and columns etc all at the same time 
  B=A(idx{:});
elseif r==-1
  np=size(ifinv,2);
  dm=sum(size(A)>1);
  B=squeeze(reshape(zeros(1,np.^dm),[1 repelem(np,dm)]));
  % commenting this out to prevent offending trilos(): B(:)=deal(NaN);
  B(idx{:})=A;
end
