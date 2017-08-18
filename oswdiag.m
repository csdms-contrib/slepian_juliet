function oswdiag(fid,fmts,lpars,thhat,thini,scl,ts,vHxs,momx,covFHh)
% OSWDIAG(fid,fmts,lpars,thhat,thini,scl,ts,vHxs,momx,covFHh)
%
% Writes an entry in the DIAGN file for the Olhede & Simons (2013) suite
%
% INPUT:
%
% fid              File identifier
% fmts             A cell array with format strings, from OSOPEN
% K                Number of independent wavenumbers, equal to 
%                  length(k(~~k))/2 [entire plane]
% lpars            lpars{1} is the log-likelihood at the estimate, from FMINUNC/FMINCON
%                  lpars{2} is the score at the estimate, from FMINUNC/FMINCON
%                  lpars{3} is the Hessian at the estimate, from FMINUNC/FMINCON
%                  lpars{4} is the exit flag, from FMINUNC/FMINCON
%                  lpars{5} is the output structure, from FMINUNC/FMINCON
% thhat,thini,scl  Estimates, initial values, scales
% ts               Optimization timing
% vHxs             Spatial (sample) variance
% momx             Moments of the quadratic portion of the likelihood
% covFHh           A single covariance matrix for the estimate, watch the calling
%                  function (Analytical Fisher-based? Analytical
%                  blurred-Hessian-based? Numerical Hessian based? Evaluated
%                  where exactly? Only one is requested as input and
%                  written to file)
%
% SEE ALSO: OSRDIAG
%
% Last modified by fjsimons-at-alum.mit.edu, 08/18/2017

% There may be one, there may be two sample variances incoming
fmtx=repmat(' %9.3e',1,length(vHxs));

% First line: The initial guess, scaled back to proper units
fprintf(fid, fmts{1},                     thini.*scl      );
% Second line: The estimate, scaled back to proper units, the sample variance
fprintf(fid,[fmts{1}(1:end-2) fmtx '\n'],[thhat.*scl vHxs]);
% Third line: Time spent, Exit flag, Iterations
% Fourth line: Likelihood, First-order optimality , moments
% Fifth line: Scales
% Remaining lines: Numerical Score, Numerical Hessian, Covariance
fprintf(fid,fmts{3},...
	round(ts),lpars{4},lpars{5}.iterations,...
	lpars{1},lpars{5}.firstorderopt,momx,...
	scl,...
        lpars{2},trilos(lpars{3}),trilos(covFHh));
% Last line: Empty (only visuals, no effect on reading)
fprintf(fid,'\n');
