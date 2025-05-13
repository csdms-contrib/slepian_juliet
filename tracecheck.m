function tracecheck(L,A,mv,tol,meth)
% TRACECHECK(L,A,mv,tol,meth)
%
% Checks the trace (sum of the eigenvalues), or the sum of the squares of
% the eigenvalues, of a matrix quadratic form L'*A*L, against an
% analytically derived expression. No output for pass, error for fail.
%
% INPUT:
%
% L       The cell array with the various first matrices
% A       The cell array with the various second matrices
% mv      This which should be equal to twice the negative trace or 
%         which should be equal to the sum of the squares of the eigenvalues
% tol     The negative of the tolerance exponent
% meth    1 Checks the sum of the eigenvalues
%         2 Checks the sum of the squares of the eigenvalues
%
% SEE ALSO:
%
% MAOSL
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/20/2018
% Last modified by olwalbert-at-princeton.edu, 02/26/2025

for j=1:length(A)
  % Checks a random wavenumber
  randi=ceil(rand*length(A{j}));

  % These are matrices out of MAOSL
  mrand=mv{j}(min(randi,length(mv{j})));
  Arand=[A{j}(randi,1) A{j}(randi,2) ; A{j}(randi,2) A{j}(randi,3)];
  % This is a lower-triangular matrix out of a Cholesky decomposition
  Lrand=[L(randi,1)    0             ; L(randi,2)    L(randi,3)];
  
  switch meth
   case 1
    % Eq. (122) in doi: 10.1093/gji/ggt056
    chekit=trace(Lrand'*Arand*Lrand);
    chekot=-2*mrand;
   case 2
    chekit=sum([eig(Lrand'*Arand*Lrand)].^2);
    chekot=mrand;
   otherwise
    error('Specify valid method')
  end

  if ~isnan(mrand) && ~isnan(chekit) && mrand~=0
    % This may by chance be at zero wavenumber in which case we skip the
    % test. We also see poor behavior here for the nu->Inf special case where mv is
    % always zero and relac will be NaN. Check RELATIVE accuracy, do not return
    % an output to DIFER!  % Maybe should only check the leading terms for the
    % tiny ones (for D) since that's where the finite precision is going to be
    % pretty bad
    relac=(chekit-chekot)/chekot;
    difer(relac,tol,[],NaN)
  else
    if j==1
      disp(sprintf('The %ist entity not being checked',j))
    elseif j==2
      disp(sprintf('The %ind entity not being checked',j))
    elseif j==3
      disp(sprintf('The %ird entity not being checked',j))
    else
      disp(sprintf('The %ith entity not being checked',j))
    end
  end
end
