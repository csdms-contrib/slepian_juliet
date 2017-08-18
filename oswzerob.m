function oswzerob(fid,th0,params,options,bounds,fmts)
% OSWZEROB(fid,th0,params,options,bounds,fmts)
% 
% Writes the (b)eginning of a THZRO diagnostic file
%
% INPUT:
%
% fid        The file id (from OSOPEN)
% th0        The parameter vector (see, e.g., SIMULOS)
% params     A structure with the known constants (see, e.g. SIMULOS)
% options    The options used by the optimization procedure
% bounds     The bounds used by the optimization procedure
% fmts       Cell array with format strings from OSOPEN
%
% SEE ALSO: 
%
% OSWZEROE, OSRZERO, OSWDIAG, OSRDIAG
%
% Last modified by fjsimons-at-alum.mit.edu, 08/18/2017

% Commit the truth to file
fprintf(fid,'%s\n','the true parameter vector');
fprintf(fid,fmts{1},th0);

% Commit the parameters of the experiment to file
fprintf(fid,'%s\n','the fixed experimental parameters');

% Rather, these need to be ordered to yield to the format
fulls={'DEL','g','z2','dydx','NyNx','blurs','kiso','quart'};
[~,i]=ismember(fulls,fieldnames(params));
jk=struct2cell(params);
fprintf(fid,fmts{2},[jk{i(~~i)}]);

% Convert the bounds to something printable
fprintf(fid,'%s\n','the bounds, if any');
if ~isempty(bounds)
  struct2str(cell2struct(bounds,...
		    {'A',  'B'  ,... % Linear Inequalities
		     'Aeq','Beq',... % Linear Equalities
		     'LB',...        % Lower Bounds
		     'UB',...        % Upper Bounds
		     'NONLCON'},...  % Nonlinear Inequalities
			 2),fid);
else
  % We need at least one colon on the next line
  struct2str(cell2struct({'None'},{'Bounds'},2),fid)
end

% Commit the parameters of the optimization to file
fprintf(fid,'%s\n','the optimization options');
struct2str(options,fid)
