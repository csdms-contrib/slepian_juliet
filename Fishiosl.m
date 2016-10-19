function F=Fishiosl(k,th)
% F=FISHIOSL(k,th)
%
% Calculates the entries in the Fisher matrix of Olhede & Simons (2013) 
% for the SINGLE-FIELD Matern model, post wavenumber averaging. No
% blurring is possible here, since no data are involved and we work from
% analytical expressions, see LOGLIOSL.
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The three-parameter vector argument [not scaled]:
%          th(1)=s2   The first Matern parameter, aka sigma^2 
%          th(2)=nu   The second Matern parameter 
%          th(3)=rho  The third Matern parameter 
%
% OUTPUT:
%
% F      The 6-column Fisher-k matrix, 'celled' in this order:
%          [1] Fs2s2   [2] Fnunu  [3] Frhorho
%          [4] Fs2nu   [5] Fs2rho [6] Fnurho [not scaled]
%
% SEE ALSO: 
%
% COVTHOSL, HESSIOSL, TRILOS, TRILOSI
% 
% EXAMPLE:
% 
% [~,th0,p,k]=simulosl([],[],1);
% F=Fishiosl(k,th0);
%
% Last modified by fjsimons-at-alum.mit.edu, 10/19/2016

defval('xver',1)

% Usually we remove the zero wavenumber from consideration
if ~isempty(k(~k))
  warning(sprintf('%s zero wavenumber',upper(mfilename))); 
end

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% First compute the "means" parameters
m=mAosl(k,th,xver);

% Then compute the various entries in the order of the paper
lk=length(k(:));

% Some of them depend on the wave vectors, some don't
F=cellnan([npp 1],[1 repmat(lk,1,5)],repmat(1,1,6));

% Fthsths, eq. (A60)
for j=1:3
  F{j}=m{j}.^2;
end

% All further combinations, eq. (A60)
jcombo=nchoosek(1:np,2);
for j=1:length(jcombo)
  F{np+j}=m{jcombo(j,1)}.*m{jcombo(j,2)};
end

% All the verification has already happened in the subroutines

% Take the expectation and put the elements in the right place
for ind=1:length(F)
   F{ind}=nanmean(F{ind});
end

