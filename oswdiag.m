function oswdiag(fid,fmts,lpars,thhat,thini,scl,ts,vHxs,covX)
% OSWDIAG(fid,fmts,lpars,thhat,thini,scl,ts,vHxs,covX)
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
%                  lpars{6} not written out here, but rather by OSWZEROB
%                  lpars{7} not written out here, but rather by OSWZEROB
%                  lpars{8} the residual moment statistics used for model testing 
% thhat,thini,scl  Estimates, initial values, scales
% ts               Optimization timing
% vHxs             Spatial (sample) variance
% covX             A certain covariance matrix estimate for the estimate
%
% SEE ALSO: OSRDIAG, OSRZERO
%
% Last modified by fjsimons-at-alum.mit.edu, 06/25/2018

% There may be one, there may be two sample variances incoming
fmtx=repmat(' %9.3e',1,length(vHxs));

% First line: The initial guess, scaled back to proper units
fprintf(fid, fmts{1},                     thini.*scl      );
% Second line: The estimate, scaled back to proper units, the sample variance
fprintf(fid,[fmts{1}(1:end-2) fmtx '\n'],[thhat.*scl vHxs]);
% Third line: Time spent, Exit flag, Iterations
% Fourth line: Likelihood, First-order optimality, moments
% Fifth line: Scales
% Remaining lines: Numerical Score, Numerical Hessian, Covariance
fprintf(fid,fmts{3},...
	round(ts),lpars{4},lpars{5}.iterations,...
	lpars{1},lpars{5}.firstorderopt,lpars{8},...
	scl,...
        lpars{2},trilos(lpars{3}),trilos(covX));
% Last line: Empty (only visuals, no effect on reading)
fprintf(fid,'\n');
