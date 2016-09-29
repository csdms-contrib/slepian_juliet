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
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% Hk       A complex matrix of Fourier-domain observations
%
% OUTPUT:
%
% gammak   A 3-column vector with the wavenumbers unwrapped, containing
%          the five partials of the likelihood function as columns.
%
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2015

% The number of parameters to solve for
np=length(th);

% First get the special matrices etc.
%[m,A]=mAosl(k,th);
m=mAosl(k,th);

% Then get the power spectrum
S=maternos(k,th);

% Now compute the score properly speaking
gammak=nan(length(k(:)),np);
for j=1:np
  % The below two are alternative formulations
  % gammak(:,j)=-m{j}-hformos(S,A{j},Hk);
  % gammak(:,j)=-m{j}+hformos(S,m{j},Hk);
  gammak(:,j)=-m{j}.*[1-abs(Hk).^2./S];
end
