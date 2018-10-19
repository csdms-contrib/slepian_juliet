function [k,dci,dcn,kx,ky]=knums(params,doit)
% [kor,dci,dcn,kx,ky]=KNUMS(params)
% [kblur,kzero,dcn,kx,ky]=KNUMS(params,1)
% 
% A partial-output interface to KNUM2 for the Olhede & Simons suite
%
% INPUT:
%
% params   A structure with AT LEAST these constants
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
% doit     1 Actually USE the params.blurs value to interpolate the k-grid 
%          (Anything goes for this parameter, only the presence counts)
% 
% OUTPUT (either THREE, or TWO, depending on the input, in addition to kx,ky):
%
% kor       The wavenumber matrix (the norm of the wave vectors)
% dci       The (m,n) indices to the DC component [at floor(dim/2)+1]
% dcn       The (m,n) indices to the components at exactly the Nyquist
%           (if any: not for odd-length dimensions; one-sided for even)
% 
% kblur     The interpolated wavenumber axis all ready to be blurred
% kzero     The running index of the zero-wavenumber in the unblurred matrix
% dcn       The (m,n) indices to the components at the Nyquist in kblur
%           (if any: not for odd-length dimensions; one-sided for even)
%
% kx,ky     The components of the wave vector
%
% Last modified by fjsimons-at-alum.mit.edu, 08/24/2017

% Extract the variables explicitly from this structure
NyNx=params.NyNx;
dydx=params.dydx;
blurs=params.blurs;

if nargin==1 || [nargin==2 && [blurs == 0 || blurs == 1]]
  % We usually run KNUM2 once focused on the output grid
  [k,kx,ky,dci,dcn]=knum2(NyNx,[(NyNx(1)-1)*dydx(1) (NyNx(2)-1)*dydx(2)]);
else
  if blurs<0
    error('You should be running BLUROSY, not BLUROS!')
  end
  
  % Proposed new dimensions
  bNyNx=blurs*NyNx;

  % Protect against parity CHANGE, which messes with BLUROS Nyquist notions
  pp=mod(NyNx,2)~=mod(bNyNx,2); 
  bNyNx=bNyNx+pp;

  % And then we run KNUM2 again to do the blurring later
  [k,kx,ky,~,dcn]=knum2(bNyNx,[(bNyNx(1)-1)*dydx(1) (bNyNx(2)-1)*dydx(2)]);
  % But then we still will want the RUNNING index of the zero in the
  % UNBLURRED matrix without running this same function again
  dci=sub2ind(NyNx,floor(NyNx(1)/2)+1,floor(NyNx(2)/2)+1);
end
