function blurcheck(Sbar,params)
% BLURCHECK(Sbar,params)
%
% Performs a check on the blurring of a spectral matrix. Checks Hermiticity
% in the univariate case, when Sbar is a column of length prod(params.NyNx),
% and the eigenvalues of Sbar for being real and positive in the bivariate
% case, in which case Sbar has dimensions prod(params.NyNx)-by-3.
%
% INPUT:
%
% Sbar    The unwrapped blurred spectral matrix being checked
%         For univariate spectra, it's a column vector. 
%         For bivariate spectra, it's [S_XX(k(:)) SXY(k(:)) SYY(k(:))]
% params  Parameters of this experiment, the ones that are needed are:
%         NyNx  number of samples in the y and x directions
%
% Last modified by fjsimons-at-alum.mit.edu, 11/02/2016

% Target dimensions, the original ones
NyNx=params.NyNx;

% Uni-variate case
if size(Sbar,2)==1
  % Check Hermiticity
  hermcheck(reshape(Sbar,NyNx));
end

% Bi-variate case
if size(Sbar,2)==3
  % Number of random tries
  ntry=10;
  for kindex=1:ntry
    % Some random wavenumber entry
    kk=max(round(rand*prod(NyNx)),1);
    % Check the eigenvalues of the little matrix at that wavenumber
    egos=eig([Sbar(kk,1) Sbar(kk,2) ; Sbar(kk,2) Sbar(kk,3)]);
    if ~all(egos>0)
      disp(sprintf('%s: some NEGATIVE eigenvalues',upper(mfilename)))
    end      
    if ~isreal(egos) 
      disp(sprintf('%s: some IMAGINARY eigenvalues',upper(mfilename)))
      cpxity(kindex)=100*mean(abs(imag(egos(:))))./mean(abs(real(egos(:))));
    end
  end
  % Attempt a brief report
  try
    disp(sprintf(...
        '%s: maximum IM/RE percentage out of %i tried is %5.0e%s',...
        upper(mfilename),ntry,max(cpxity),'%'))
  catch
    %disp(sprintf('%s: no problems found',upper(mfilename)))
  end
end
