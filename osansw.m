function [answ,answs]=osansw(th0,covX,E,v)
% [answ,answs]=OSANSW(th0,covX,E,v)
%
% Constructs some reporting strings for MLEOS, MLEROS, MLEROS, MLEOSL, 
% thus for correlated and uncorrelated loading scenarios as well as the
% single-field estimate (which does not require the third and fourth input.)
%
% INPUT:
%
% th0        True parameter vector that you want quoted
% covX       The estimation covariance that you want quoted
% E          Young's modulus [Pa] (for dual fields, using DTOTE)
% v          Poisson's ratio ] (for dual fields, using DTOTE)
%
% OUTPUT:
%
% answ       A cell with the "answers" if you will
% answs      Appropriate formatting strings for a plot title
%
% SEE ALSO: OSDISP, DTOTE
%
% Last modified by fjsimons-at-alum.mit.edu, 10/17/2017

defval('th0',[])
defval('covX',[])
defval('E',[])
defval('v',[])

if all(isempty([th0 covX(:)' E v]))
  answ={''};
  answs='%s';
end

% The plus-minus sign
pm=str2mat(177);

% Dual-field analysis with no business for E or v
if ~isempty(E) && ~isempty(v)
  % Conversion to Te [m]
  Te=DtoTe(th0(1),E,v);
  
  % This will come awfully close. We remember to put in the mean!
  stdTe=std(DtoTe(th0(1)+randn(10000,1)*sqrt(covX(1,1)),E,v));
  stdf2=sqrt(covX(2,2));
  
  % Here is the "delta method", which is approximate
  stdTedelta=sqrt(covX(1,1))/th0(1)^(2/3)/3*[12*(1-v^2)/E]^(1/3);

  disp(sprintf('Stdv of the Te predicted by simulation:  %3.2f km',  stdTe/1000))
  disp(sprintf('Stdv of the Te predicted by delta meth:  %3.2f km\n',stdTedelta/1000))
  
  % True parameters and their purported uncertainties
  answ{1}=sprintf('E = %g',E);
  answ{2}=sprintf('%s = %g','\nu',v);
  answ{3}=sprintf('T_e = %4.1f %s %3.1f km',...
		  round(Te/1e3*10)/10,pm,round(stdTe/1e3*10)/10);
  answ{4}=sprintf('f^2 = %g %s %g',...
		  th0(2),pm,round(stdf2*1e3)/1e3);
end
  
% Some common strings
s2s='%s^2 = %6.0f %s %6.0f';
nus='%s = %4.2f %s %4.2f';
rhs='%s = %i %s %i';

if length(th0)==3
  % Single-field analysis (MLEOSL etc)
  stds2  =sqrt(covX(1,1));
  stdnu  =sqrt(covX(2,2));
  stdrho =sqrt(covX(3,3));

  answ{1}=sprintf(s2s,'\sigma',th0(1),pm,stds2);
  answ{2}=sprintf(nus,'\nu',   th0(2),pm,stdnu);
  answ{3}=sprintf(rhs,'\rho',  round(th0(3)),pm,round(stdrho));

  answs=['%s ; %s ; %s'];
elseif length(th0)==5
  % Uncorrelated dual-field analysis (MLEOS etc)
  stds2  =sqrt(covX(3,3));
  stdnu  =sqrt(covX(4,4));
  stdrho =sqrt(covX(5,5));
  
  answ{5}=sprintf(s2s,'\sigma',th0(3),pm,stds2);
  answ{6}=sprintf(nus,'\nu',   th0(4),pm,stdnu);
  answ{7}=sprintf(rhs,'\rho',  round(th0(5)),pm,round(stdrho));

  answs=['%s ; %s ; %s ; %s\n %s ; %s ; %s'];
elseif length(th0)==6
  % Correlated dual-field analysis (MLEROS etc)
  stdr   =sqrt(covX(3,3));
  stds2  =sqrt(covX(4,4));
  stdnu  =sqrt(covX(5,5));
  stdrho =sqrt(covX(6,6));
  
  answ{5}=sprintf('r = %g %s %g',th0(3),pm,round(stdr*1e3)/1e3);
  answ{6}=sprintf(s2s,'\sigma',	 th0(4),pm,stds2);
  answ{7}=sprintf(nus,'\nu',	 th0(5),pm,stdnu);
  answ{8}=sprintf(rhs,'\rho',	 round(th0(6)),pm,round(stdrho));
  
  answs=['%s ; %s ; %s ; %s ; %s\n %s ; %s ; %s'];
end
