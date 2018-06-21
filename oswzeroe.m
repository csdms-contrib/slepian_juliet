function oswzeroe(fid,scl,avhsz,good,F0,covF0,fmti)
% OSWZEROE(fid,scl,avhsz,good,F0,covF0,fmti)
% 
% Writes the (e)nd of a THZRO diagnostic file
%
% INPUT:
%
% fid          The file id (from FOPEN)
% scl          The scaling factors for the Fisher matrix 
% avhsz        The full-form sample-average scaled Hessian matrix
% good         The number of samples over which this was averaged
% F0           The full-form scaled Fisher matrix, at the truth
% covF0        The full-form theoretical covariance matrix, based on F0
% fmti         Strings containing formatting instructions    
%
% SEE ALSO: 
%
% OSWZEROB, OSRZERO, OSWDIAG, OSRDIAG
%
% Last modified by fjsimons-at-alum.mit.edu, 08/21/2017

% Print the scaling of the theoretical values
fprintf(fid,'%s\n','The scaling factors');
fprintf(fid,fmti{4},scl);

% Now print the unscaled theoretical covariance to file also
fprintf(fid,'%s\n','The covariance from the Fisher matrix at the truth');
fprintf(fid,fmti{1},trilos(covF0));

% Print the scaled Fisher matrix, the expected value of the Hessian
fprintf(fid,'%s\n','The unblurred scaled Fisher matrix at the truth');
fprintf(fid,fmti{2},trilos(F0));

% Print the observed scaled average of the Hessians (over last set of runs)
fprintf(fid,'%s\n',...
        sprintf('The scaled numerical Hessian matrix averaged over whichever %i runs last output',...
                good));
fprintf(fid,fmti{2},trilos(avhsz));
