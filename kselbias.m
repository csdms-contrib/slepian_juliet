function kselk=kselbias(k,p)
% kselk=KSELBIAS(k,p)
%
% Finds the ratio between the number of selected wavenumbers and the
% available wavenumbers to enable efficiency calculations down the line.
%
% INPUT:
%
% k        Original wavenumbers
% params   A parameter structure with AT LEAST this constant:
%            kiso   wavenumber beyond which we are NOT considering
%
% OUTPUT:
%
% kselk    The ratio of the number of retained over the available entries
%
% SEE ALSO:
%
% LOGLIOSL, VARBIAS
%
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2018

defstruct('p','kiso',NaN)

if any(~isnan(p.kiso))
  kselk=sum(k(:)<=p.kiso)/prod(size(k));
else
  kselk=1;
end




