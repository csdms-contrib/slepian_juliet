function [k,dci,dcn,kx,ky]=knums(params,doit)
% [kor,dci,dcn,kx,ky]=KNUMS(params)
% [kblur,kzero,dcn,kx,ky]=KNUMS(params,1)
% 
% Produces a grid of wavenumbers suitable for the spectral analysis
% of spatial data as part of the Olhede & Simons suite
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
% OUTPUT (either THREE, or TWO, depending on the input, IN ADDITION to kx,ky):
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
% SEE ALSO:
%
% KNUM2, which is called by this function, BLUROS
%
% Last modified by fjsimons-at-alum.mit.edu, 03/02/2022

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

  % And then we run KNUM2 again to do the blurring later
  [k,kx,ky,~,dcn]=knum2(bNyNx,[(bNyNx(1)-1)*dydx(1) (bNyNx(2)-1)*dydx(2)]);
  % But then we still will want the RUNNING index of the zero in the
  % UNBLURRED matrix without running this same function again, so
  % produce that under the fake name that gets used on the inside only
  dci=sub2ind(NyNx,floor(NyNx(1)/2)+1,floor(NyNx(2)/2)+1);
end

% SEE BLUROSY TO FIGURE OUT THE REVERSE SUBSAMPLING!

% CONSIDERATIONS THAT WERE AT ONE TIME PART OF TESTING
%
% We should indeed be able to use ANY refinement grid. 
% Should we protect against odd->even parity CHANGE? 
% pp=mod(NyNx,2)~=mod(bNyNx,2);
% As gleaned from BLUROS('demo1') there remains a slight preference for
% combining odd dimensions with odd refinements, it's just got to do with
% the way they funky Fejk works out. It's not an interpolation issues -
% we're subsampling everything now. The trick is to keep the product of
% the grid size and the refined factor as high as possible .
% e o -> e v no parity change of original
% o o -> o v no parity change of original
% e e -> e v no parity change of original
% o e -> e x no parity change of original
% This last case gave us so much of the trouble that helped us track down
% all the other mistakes. Not it's just slightly less accurate, that's all. 
