function varargout=randgpn1(k,dc,dcn,xver)
% [Zk,Zx]=randgpn1(k,dc,dcn,xver)
% 
% Returns a set of (complex proper) normal variables suitable for
% IFFT, most notably this works for even and/or odd length vectors
%
% INPUT:
%
% k       A wavenumber vector (DC component in center)
% dc      The index to the DC component in k
% dcn     The index to the Myquist components in k
%         --> Note that these three are straight out of KNUM1, and that
%         the actual wavenumbers are not used, only size(k) is needed. 
% xver    1 Checks the Hermiticity of the result by inverse transformation
%         0 No such check is being conducted 
%
% OUTPUT:
%
% Zk      A size(k) vector with (complex proper) normal Fourier variables
% Zx      A size(k) vector with a real-valued random field
%
% EXAMPLE 1
%
% [k,kx,dci,dcn]=knum1(round(rand*100),100);
% [Zk,Zx]=randgpn1(k,dci,dcn); plot(Zx)
%
% EXAMPLE 2 
% 
% [k,kx,dci,dcn]=knum1(200,100);
% F=abs(randn*10); [Zk,Zx]=randgpn1(k,dci,dcn); Zf=ifft(ifftshift(Zk.*exp(-F*k)));
% plot(Zf); title(sprintf('exp(%3.3gk)',-F));
%
% SEE ALSO: KNUM2, KNUM1
%
% Last modified by fjsimons-at-alum.mit.edu, 03/18/2019

defval('xver',0)

% Make a receptacle with the dimension of the wavenumbers
n=length(k);
Zk=zeros(1,n); clear k

% Determine the parity
nodd=mod(n,2);

% Define the ranges to pick out the halves
lhn=1:dc;
rhn=dc+1:n;

% Now make some complex proper normal random variables of half variance
ReZk=randn(1,dc)/sqrt(2);
ImZk=randn(1,dc)/sqrt(2);
% And make a handful real normal random variables of unit variance
realZk=randn(2,1);

% Fill the entire left half plane with these random numbers
lh=ReZk+sqrt(-1)*ImZk;
Zk(1,lhn)=lh;

% Fill the center with a real, where the DC component goes 
Zk(1,dc)=realZk(1);

% Fill the Nyquist with a real, if we capture them exactly 
for ind=1:size(dcn,1)
  Zk(dcn)=realZk(1+ind);
end

% And now symmetrize the right half plane so that output is Hermitian
Zk(1,rhn)=conj(fliplr(lh(2-nodd:dc-1)));

if nargout>1 || xver==1
  Zx=ifft(ifftshift(Zk));
  % You could now check that the IFFT2 is real (don't forget the "2"!!)
  if ~isreal(Zx) ; error(sprintf('Not Hermitian by %5g',mean(imag(Zx(:))))); end
  % You can perform this check more directly as well
  hermcheck(Zk)
else
  Zx=NaN;
end

% Output
varns={Zk,Zx};
varargout=varns(1:nargout);
