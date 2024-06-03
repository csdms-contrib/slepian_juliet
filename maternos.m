function varargout=maternos(k,th,dth,d,meth)
% SkordSkdth=MATERNOS(k,th,dth,d,meth)
%
% Calculates the three-parameter isotropic d-dimensional Matern
% spectral density, which is dependent on d, also the counterpart to the
% isotropic correlation. See Olhede & Simons (2013), doi: 10.1093/gji/ggt056.x 
% Additionally can calculate the partial derivatives of the spectral density
% with respect to the Matern parameters.
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m]
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2   The first Matern parameter [variance in units^2]
%          nu   The second Matern parameter [differentiability]
%          rho  The third Matern parameter [range in m]
% d        The dimensionality [default is 2]
% meth     1 straight up calculation
%          2 analytical simplification for special values of th(end-1)
% dth      When empty, the spectral density is calculated; otherwise,
%          1, 2, or 3 specifies which element of th gets differentiated
%
% OUTPUT:
%
% SkordSkdth  A column vector with all the wavenumbers unwrapped, possibly
%             calculated using a simplification of nu, OR returns the DERIVATIVE in
%             the dth element of th; simplifications and derivatives are included 
%             to mirror MATERNOSY, however, these do not feed into another program
%
% SEE ALSO:
%
% MATERNPRC, MATERNOSY, MATERNOSP
%
% EXAMPLE:
%
% [~,th0,p]=simulosl; difer(maternos(0,th0)/th0(1)-pi*th0(3)^2/4,7,[],NaN)
% 
% Periodogram for the inifinite smoothness, Gaussian, squared exponential case:
%
% th0(2)=Inf; k=knums(p); k=k(:); Sk=maternos(k,th0);
% figure(); imagesc(log10(v2s(Sk,p))); axis image
%
% Confirmation that the special case simplifications are equivalent to the
% straight calculation:
%
% nus=[1/3,1,1/2,3/2,5/2,7/2,15/2];th=[1 nus(randi(length(nus),1)) 1]
% Sk=maternos(k,th,[],2,1); Sk2=maternos(k,th,[],2,2); all(Sk==Sk2)
%
% Calculate and visualize the partial derivatives: 
%
% dth=randi(3,1); dSkdth=maternos(k,th,[],2,dth); 
% figure(); imagesc(log10(v2s(dSkdth,p))); axis image
%
% Last modified by fjsimons-at-alum.mit.edu, 06/03/2024
% Last modified by olwalbert-at-princeton.edu, 06/03/2024

% Defaults
defval('d',2)
defval('meth',1)
defval('dth',[])

% The Matern parameters are always the last three elements of the TH input 
s2=th(end-2);
nu=th(end-1);
rh=th(end  );

% Adjust for dimensionality by specialization
switch d
 case 2
  pd=nu;
 otherwise
  pd=gamma(nu+d/2)/gamma(nu);
end

% Calculate the spectral density?
if isempty(dth)
    % Check whether we are calculating the spectral density for the case that
    % nu approaches infinity
    if isinf(nu)
        % Squared exponential
        % Calculate the d-dimensional squared exponential spectral density. 
        % This analytic form of the Matern spectral density was solved for 
        % using Eq. 6.1.46 of Abramowitz & Stegun (1965) and L'Hopital's rule
        % for indeterminate limits.
        Sk=s2*(pi*rh^2/4)^(d/2)*exp(-pi^2*rh^2*k(:).^2/4); 
    else
        % Calculate the denominator in the spectral density
        avark=4*nu/pi^2/rh^2+k(:).^2;
        % Switch the calculation method
        switch meth 
            case 1
              % The full evaluation, calculate the d-dimensional spectral
              % density
              Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2);
            case 2
              % By selecting the second method of calculation, we are seeking
              % to evaluate Sk from a simplified analytic expression of the
              % d-dimensional Matern spectral density for a special value of
              % nu
              if nu==1/3
                  % von Karman
                  if d==2
                    Sk=s2*2^(2/3)*pi*rh^2/(4+3*(pi*rh*k(:)).^2).^(4/3); 
                  else
                    Sk=gamma(1/3+d/2)/gamma(1/3)*s2*...
                        (4/(4+3*(pi*rh*k(:)).^2)).^(1/3)*...
                        (3*pi*rh^2/(4+3*(pi*rh*k(:)).^2))^(d/2);
                  end
              elseif nu==1/2
                  % Exponential
                  if d==2
                      Sk=s2*pi*rh^2/(sqrt(2)*(2+(pi*rh*k(:)).^2).^(3/2));
                  else
                      Sk=gamma((d+1)/2)*s2*(2/(pi^(1-d)))^(1/2)*rh^d*...
                          (2+(pi*rh*k(:)).^2).^(-(d+1)/2);
                  end
              elseif nu==1
                  % Whittle
                  if d==2
                    Sk=s2*4*pi*rh^2/(4+(pi*rh*k(:)).^2).^2; 
                  else
                    Sk=gamma(1+d/2)*s2/pi^(d/2)*(2/(pi*rh))^2*...
                        ((pi*rh)^2/(4+(pi*rh*k(:)).^2)).^(1+d/2);
                  end
              elseif nu==3/2
                  % Second-order autoregressive
                  if d==2
                    Sk=9*sqrt(6)*s2*pi*rh^2/(6+(pi*rh*k(:)).^2).^(5/2);
                  else
                    Sk=gamma((d+3)/2)*12*(6*pi)^(1/2)*s2*rh^2/...
                        (6+(pi*rh*k(:)).^2).^(5/2);
                  end
              elseif nu==5/2
                  % Third-order autoregressive
                  if d==2
                    Sk=250*sqrt(10)*s2*pi*rh^2/(10+(pi*rh*k(:)).^2).^(7/2);
                  else
                    Sk=gamma((d+5)/2)*s2*400/3*(10)^(1/2)*pi^(d/2-1/2)*rh^d*...
                        (10+(pi*rh*k(:)).^2).^(-(d+5)/2);
                  end
              elseif mod(nu-1/2,1)==0
                  % A general half-integer case beyond 1/2, 3/2, or 5/2,
                  % nu=n+1/2 for n = 1, 2, 3, ...
                  n=nu-1/2;
                  if d==2
                    Sk=s2*pi*rh^2/4*((4*n+2)/(4*n+2+(pi*rh*k(:)).^2)).^(n+3/2);
                  else
                    Sk=gamma(n+(d+1)/2)/gamma(n+1/2)*s2*((4*n+2)/...
                        (4*n+2+(pi*rh*k(:)).^2)).^(n+1/2)*...
                        (pi*rh^2/(4*n+2+(pi*rh*k(:)).^2)).^(d/2);
                  end
              else
                  % However, if the nu provided to MATERNOS is not one of the
                  % special values that we have considered so far, we should
                  % throw an error
                  error('This is not a special case of nu. Ask for meth=1.')
              end 
        end 
    end
% Or, calculate a partial derivative?
else
   % Calculate the partial derivative wrt the requested parameter index
   % disp(sprintf('Calculating %i%s parameter derivative',dth,ith(dth)))
   % Abuse of nomenclature, now you're getting a specific derivative
   if dth==1
       % The partial derivative of Sk with respect to the variance, dSkds2
       Sk=pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2);
   elseif dth==2
       % The partial derivative of Sk with respect to the smoothness, dSkdnu; 
       % see Eq. 6.3.1 of Abramowitz & Stegun (1965) for the derivative of the 
       % Gamma function and 
       Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2)*...
          (psi(nu+d/2)-psi(nu)+log(nu*4/(pi*rh)^(2))-log(avark)-...
          (4*nu+2*d)/(4*nu+(pi*rh*k(:)).^2)+1);
   elseif dth==3
       % The partial derivative of Sk with respect to the range, dSkdrh
       Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2-1)*...
          (2*nu/rh)*(2*d/(pi*rh)^2-k(:).^2);
   else
       error('Not a valid partial derivative index. Ask for dth=1,2,or 3.') 
   end
end
% Serve output
varns={Sk};
varargout=varns(1:nargout);
