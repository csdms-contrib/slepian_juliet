function tracecheck(L,A,mv,tol,meth)
% TRACECHECK(L,A,mv,tol,meth)
%
% Checks the trace (sum of the eigenvalues), or the sum of the squares of
% the eigenvalues, of a matrix L'*A*L, against an analytically derived
% expression.
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
% If it passes, no output, if it doesn't, get an error
% 
% Last modified by fjsimons-at-alum.mit.edu, 10/08/2014

for j=1:length(A)
  % Checks a random wavenumber
  randi=ceil(rand*length(A{j}));

  mrand=mv{j}(min(randi,length(mv{j})));
  Arand=[A{j}(randi,1) A{j}(randi,2) ; A{j}(randi,2) A{j}(randi,3)];
  Lrand=[L(randi,1)    0       ; L(randi,2) L(randi,3)];
  
  switch meth
   case 1
    chekit=trace(Lrand'*Arand*Lrand);
    chekot=-2*mrand;
   case 2
    chekit=sum([eig(Lrand'*Arand*Lrand)].^2);
    chekot=mrand;
   otherwise
    error('Specify valid method')
  end
  
  if ~isnan(mrand) && ~isnan(chekit)
    % This may by chance be at zero wavenumber in which case we skip the
    % test. Check RELATIVE accuracy, do not return an output to difer!
    % Maybe should only check the leading terms for the tiny ones (D)
    % since that's where the finite precision is going to be pretty bad
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
