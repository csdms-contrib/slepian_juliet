function oswzeroe(fid,scl,avH,good,F,covF,fmti)
% OSWZEROE(fid,scl,avH,good,F,covF,fmti)
% 
% Writes the (e)nd of a THZRO diagnostic file
%
% INPUT:
%
% fid          The file id (from FOPEN)
% scl          The scaling factors for the Fisher matrix 
% avH          The full-form sample-average scaled Hessian matrix
% good         The number of samples over which this was averaged
% F            The full-form scaled Fisher matrix, at the truth
% covF         The full-form theoretical covariance matrix, based on F
% fmti         Strings containing formatting instructions    
%
% SEE ALSO: 
%
% OSWZEROB, OSRZERO, OSWDIAG, DIAGNOS
%
% Last modified by fjsimons-at-alum.mit.edu, 11/15/2016

% Print the scaling of the theoretical values
fprintf(fid,'%s\n','The scaling factors');
fprintf(fid,fmti{4},scl);

% Now print the unscaled theoretical covariance to file also
fprintf(fid,'%s\n','The covariance from the Fisher matrix at the truth');
fprintf(fid,fmti{1},trilos(covF));

% Print the scaled Fisher matrix, the expected value of the Hessian
fprintf(fid,'%s\n','The unblurred scaled Fisher matrix at the truth');
fprintf(fid,fmti{2},trilos(F));

% Print the observed scaled average of the Hessians (over last set of runs)
fprintf(fid,'%s\n',...
        sprintf('The scaled numerical Hessian matrix averaged over whichever %i runs last output',...
                good));
fprintf(fid,fmti{2},trilos(avH));
