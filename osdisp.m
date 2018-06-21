function varargout=osdisp(th0,thhats,nl,avhs,F0,covavhs)
% OSDISP(th0,thhats,nl,avhs,F0,covavhs)
% OSDISP(th0,params)
% [str0,str1,str2,str3]=OSDISP(...)
%
% Displays some messages to the command screen for MLEOS, MLEROS, MLEROS0
% and also SIMULOS, SIMULROS, SIMULROS0
%
% INPUT:
%
% th0        True parameter vector
% thhats     Estimated parameter vector, OR
%  params    A structure with the fixed parameter settings, AND
% nl         Number of experiments over which the average Hessian is reported
% avhs       Median numerical Hessian matrix at the estimates (after TRILOS)
% F0         Fisher matrix evaluated at at the truth (after TRILOS)
% covavhs    The covariance based on the average median Hessian (full form)
%
% OUTPUT:
%
% The strings used 
%
% Last modified by fjsimons-at-alum.mit.edu, 06/20/2018

% The necessary strings for formatting
str0='%18s';
str1='%13.5g ';
str2='%13i ';
str3='%13s ';

% Replicate to the size of the parameter vector
str1s=repmat(str1,size(th0));
str2s=repmat(str2,size(th0));
str3s=repmat(str3,size(th0));

% Don't use STRUC2ARRAY since we want them in our own order
% But see the reordering solution in OSWZEROB
if isstruct(thhats) && nargin==2
  params=thhats;
  if length(th0)>3
    disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3,1,10)),...
		 ' ','D1','D2','g','z2','dy','dx','Ny','Nx','blurs','quart'))
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,1,10)),...
                 'Parameters',params.DEL,params.g,params.z2,...
                 params.dydx,params.NyNx,params.blurs,params.quart))
  else
    disp(sprintf(sprintf('\n%s   %s ',str0,repmat(str3,1,6)),...
		 ' ','dy','dx','Ny','Nx','blurs','quart'))
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,1,6)),...
                 'Parameters',params.dydx,params.NyNx,params.blurs,params.quart))
  end
else
  % Estimated values
  disp(sprintf(sprintf('%s : %s \n',str0,str1s),...
	       'Average thhat',mean(thhats,1)))
  
  if nl==1
    % Average numerical Hessian near the estimate and Fisher matrix at the truth
    disp(sprintf(['Over %i simulation, the numerical Hessian near the estimate\n' ...
		  'and the unblurred Fisher matrix at the truth are |%4.2f|%s apart\n'],nl,...
		 1/100*round(100*mean(abs([avhs-F0]'./F0'*100))),'%'))
  else
    % Average numerical Hessian near the estimate and Fisher matrix at the truth
    disp(sprintf(['Over %i simulations, the median numerical Hessian near the estimate\n' ...
		  'and the unblurred Fisher matrix at the truth are |%4.2f|%s apart\n'],nl,...
		 1/100*round(100*mean(abs([avhs-F0]'./F0'*100))),'%'))
  end

  % Covariance, relative, empirical, and theoretical
  if size(thhats,1)>4
    disp(sprintf(sprintf('%s : %s',str0,str1s),...
		 'Observed standard deviation',std(thhats)))
  end
  disp(sprintf(sprintf('%s : %s',str0,str1s),...
	       'AvNumHes standard deviation',sqrt(diag(covavhs))))
  if size(thhats,1)>4
    disp(sprintf(sprintf('%s : %s',str0,str2s),...
		 'Perct of obs to pred stnddv',...
		 round(100*std(thhats)./sqrt(diag(covavhs)'))))
    disp(sprintf(sprintf('%s : %s\n',str0,str2s),...
		 'Observed percent stand devn',...
		 round(100*std(thhats)./th0)))
  end
end

if length(th0)==6
  disp(sprintf(sprintf('\n%s   %s ',str0,str3s),...
	       ' ','D','f2','r','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str1s),'True theta',th0))
elseif length(th0)==5
  disp(sprintf(sprintf('\n%s   %s ',str0,str3s),...
	       ' ','D','f2','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str1s),'True theta',th0))
elseif length(th0)==3
  disp(sprintf(sprintf('\n%s   %s ',str0,str3s),...
	       ' ','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str1s),'True theta',th0))
end

% Optional output
varns={str0,str1,str2,str3};
varargout=varns(1:nargout);
