function oswzerob(fid,th0,params,lpars,fmts)
% OSWZEROB(fid,th0,params,lpars,fmts)
% 
% Writes the (b)eginning of a THZRO diagnostic file
%
% INPUT:
%
% fid        The file id (from OSOPEN)
% th0        The parameter vector (see, e.g., SIMULOS)
% params     A structure with the known constants (see, e.g. SIMULOS)
% lpars      lpars{1} not written out here, but rather by OSWDIAG
%            lpars{2} not written out here, but rather by OSWDIAG
%            lpars{3} not written out here, but rather by OSWDIAG
%            lpars{4} not written out here, but rather by OSWDIAG
%            lpars{5} not written out here, but rather by OSWDIAG
%            lpars{6} the options used by the FMINUNC/FMINCON procedure
%            lpars{7} any bounds used by the  FMINUNC/FMINCON procedure
%            lpars{8} not written out here, but rather by OSWDIAG
% fmts       Cell array with format strings from OSOPEN
%
% SEE ALSO: 
%
% OSWZEROE, OSRZERO, OSWDIAG, OSRDIAG
%
% Last modified by fjsimons-at-alum.mit.edu, 05/19/2025

% Commit the truth to file
fprintf(fid,'%s\n','the true parameter vector');
fprintf(fid,fmts{1},th0);

% Commit the parameters of the SIMULATION to file
fprintf(fid,'%s\n','the fixed simulation parameters');

% Rather, these need to be ordered like this to yield to the format
fulls={'DEL','g','z2','dydx','NyNx','blurs','kiso','quart'};
[~,i]=ismember(fulls,fieldnames(params));
jk=struct2cell(params);
fprintf(fid,fmts{2},[jk{i(~~i)}]);

% Commit the parameters of the INVERSION to file
fprintf(fid,'%s\n','the fixed inversion parameters');

% blurs nugget ifinv maskhash 

% Convert the bounds to something printable
fprintf(fid,'%s\n','the bounds, if any');
if ~isempty(lpars{7})
  struct2str(cell2struct(lpars{7},...
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
struct2str(lpars{6},fid)
