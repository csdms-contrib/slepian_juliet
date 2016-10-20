function gam=gammiosl(k,th,params,Hk,xver)
% gam=GAMMIOSL(k,th,params,Hk,xver)
%
% Calculates the entries in the score matrix of Olhede & Simons (2013) 
% for the SINGLE-FIELD Matern model, post wavenumber averaging. When
% blurred, no consideration is given to the zero wavenumber, see LOGLIOSL. 
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter, aka sigma^2 
%          th(2)=nu   The second Matern parameter 
%          th(3)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%               -1 Blur likelihood using the exact procedure
%          NOTE: It's not going to be a great derivative unless you could
%          change MAOSL also. Still, the order of magnitude will be OK.
% Hk       A complex matrix of Fourier-domain observations
% xver     Excessive verification
%
% OUTPUT:
%
% gam      The derivative of the log-likelihood
%
% EXAMPLE:
% 
% p.quart=0; p.blurs=0; p.kiso=NaN; clc; [~,th0,p,k,Hk]=simulosl([],p,1);
% F=Fishiosl(k,th0); 
% G=gammiosl(k,th0,p,Hk);
% H=Hessiosl(k,th0,p,Hk);
% % On average, F and H should be close
% 
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2016

defval('xver',1)

% The number of parameters to solve for
np=length(th);
% The number of wavenumbers
lk=length(k(:));

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

% Initialize the wavenumber-dependent score vector
gammak=nan(lk,np);

% First get the special matrices etc.
if xver==0
  mth=mAosl(k,th);
elseif xver==1
  [mth,~,A]=mAosl(k,th);
end

% We need the power spectrum and its ratio to the observations 
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

% The average of Xk needs to be close to one as will be tested 
Xk=abs(Hk).^2./S;

% Now compute the score in each of the parameters
for j=1:np
  % Eq. (A53) in doi: 10.1093/gji/ggt056
  gammak(:,j)=-mth{j}.*[1-Xk];
  if xver==1
    % The below two are alternative formulations
    gkalt1(:,j)=-mth{j}-hformos(S,A{j},Hk);
    gkalt2(:,j)=-mth{j}+hformos(S,mth{j},Hk);
    % Check for numerical indistinguishability
    diferm(gkalt1(:,j),gammak(:,j));
    diferm(gkalt2(:,j),gammak(:,j));
  end
end

% Now the wavenumber averaging
gam=-nanmean(gammak);
