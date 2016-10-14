function gammak=gammakosl(k,th,params,Hk)
% gammak=GAMMAKOSL(k,th,params,Hk)
%
% Computes the first partial derivatives of the likelihood function for
% the SINGLE-FIELD isotropic Matern model in Olhede & Simons (2013).
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter, aka sigma^2 
%          th(2)=nu   The second Matern parameter 
%          th(3)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
%          NOTE: It's not going to be a great derivative unless you
%          change MAOSL also. Still, the order of magnitude will be OK.
% Hk       A complex matrix of Fourier-domain observations
%
% OUTPUT:
%
% gammak   A 3-column vector with the wavenumbers unwrapped, containing
%          the five partials of the likelihood function as columns.
%
% SEE ALSO:
%
% LKOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 10/14/2016

% The number of parameters to solve for
np=length(th);

% Initialize
gammak=nan(length(k(:)),np);

% First get the special matrices etc.
%[m,A]=mAosl(k,th);
m=mAosl(k,th);

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

% See LKOSL for the detailed explanations of these procedures
switch blurs
 case {0,1}
    S=maternos(k,th);
 otherwise
  if blurs>1
    S=bluros(maternos(knums(params,1),th),params);
  else
    S=blurosy(th,params);
  end
end

% Now compute the score properly speaking
for j=1:np
  % The below two are alternative formulations
  % gammak(:,j)=-m{j}-hformos(S,A{j},Hk);
  % gammak(:,j)=-m{j}+hformos(S,m{j},Hk);
  gammak(:,j)=-m{j}.*[1-abs(Hk).^2./S];
end

