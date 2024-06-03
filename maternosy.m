function varargout=maternosy(y,th,dth,meth)
% CyordCydth=MATERNOSY(y,th,dth,meth)
%
% Calculates d-dimensional isotropic Matern correlation, which is
% independent of d, also the counterpart to the spectral covariance.
% See Olhede & Simons (2013), doi: 10.1093/gji/ggt056.x, eq. (72)
% Additionally can calculate the partial derivatives of the isotropic Matern
% correlation with respect to the Matern parameters.
%
% INPUT:
%
% y        Lag parameter, the distance between spatial positions,
%          e.g. from XXPDIST or SSPDIST [m], can be any dimension
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2    The first Matern parameter [variance in units^2]
%          nu    The second Matern parameter [differentiability]
%          rh    The third Matern parameter [range in units]
% dth      When empty, the spatial Matern covariance is calculated; otherwise, 
%          1, 2, or 3 specifies which element of th gets differentiated
% meth     1 straight up calculation
%          2 analytical simplification for special values of th(end-1)
%
% OUTPUT:
%
% Cy       The spatial Matern covariance at all the requested lags, possibly
%          calculated using the Abramowitz & Stegun (1965) simplifications to
%          the Bessel function for half-integer values of nu, the infinite limit
%          of nu, OR returns the DERIVATIVE in the dth element of th, to feed 
%          into BLUROSY
%
% SEE ALSO:
%
% MATERNOS, MATERNOSP, MATERNPRC
%
% EXAMPLE:
%
% Discrete approximation of the integral compared to exact result:
%
% [Hx,th0,p]=simulosl; bigN=1000;
% y=linspace(0,sqrt(prod(p.dydx))*sqrt(prod(p.NyNx)),bigN);
% % S(0)=1/(2*pi)\int C(r)rdr
% [sum(y.*maternosy(y,th0))*(y(2)-y(1))/(2*pi) maternos(0,th0)]
%
% Confirmation that the special case simplifications are equivalent to the
% straight calculation:
% 
% nu=[1/3 1/2 1 3/2 5/2 31/2]; th0(2)=nu(randi(length(nu))); 
% Cy=maternosy(y,th0,[],1); Cy2=maternosy(y,th0,[],2); diferm(Cy/(Cy(1)),Cy2/Cy2(1))
%
% Calculate the infinite smoothness, Gaussian, squared exponential case:
%
% th0(2)=Inf; Cy=maternosy(y,th0);
%
% Calculate and visualize the partial derivatives: 
%
% maternosy('demo2')
%
% Demo1 provides a visual comparison of the general output of MATERNOS and
% MATERNOSY as Fourier pairs:
%
% maternosy('demo1')
%
% Last modified by fjsimons-at-alum.mit.edu, 06/03/2024
% Last modified by olwalbert-at-princeton.edu, 06/03/2024

if ~isstr(y)
  % Defaults
  defval('meth',1)
  defval('dth',[])
  % The Matern parameters are always the last three elements of the TH input 
  s2=th(end-2);
  nu=th(end-1);
  rh=th(end  );
  % The argument, make sure it is a distance
  argu=2*sqrt(nu)/pi/rh*abs(y);
  % Calculate the spatial covariance?
  if isempty(dth)
    % Check whether we are asking for the case that nu approaches
    % infinity; if so, by the regular method Cy would have Inf/NaN entries
    if isinf(nu)
        % Squared exponential 
        % This analytic form of the isotropic Matern covariance 
        % was solved from the inverse Fourier transform of the limit as
        % nu approaches infinity of the spectral density (see MATERNOS
        % for details) with use of Eq. 3.323.2 of Gradshteyn & Ryzhik (1980).
        %%%OLW: (Consistency with BLUROSY demo?)
        Cy=s2/(pi*rh)*exp(-abs(y).^2/(pi^2*rh^2));
    else
        % Switch the calculation method
        switch meth
          case 1
            % The evaluation
            Cy=2^(1-nu)*s2/gamma(nu)*argu.^nu.*besselk(nu,argu);
            % Supply the smallest arguments
            Cy(y==0)=s2;
          case 2
            % By selecting the second method of calculation, we are seeking to 
            % evaluate Cy from the simplified analytic expression of the 
            % isotropic Matern covariance for a special value of nu. The 
            % following analytic forms of half-integer nu are calculated from 
            % substitution of Eqs 10.1.9 and 10.2.15 of Abramowitz & Stegun 
            % (1965) for modified Bessel functions of half-integer orders.
            % See also 10.1093/biomet/93.4.989 for comparison.
            if nu==1/3
                % von Karman 
                Cy=s2*2^(2/3)/gamma(nu).*(2*sqrt(3)/(3*pi*rh)*abs(y)).^...
                   (nu).*besselk(nu,2*sqrt(3)/(3*pi*rh)*abs(y));
                % Compute the value at zero lag
                Cy(y==0)=s2;
            elseif nu==1/2
                % Exponential
                Cy=s2*exp(-sqrt(2)/(pi*rh)*abs(y));
            elseif nu==1
                % Whittle
                Cy=s2*2/(pi*rh)*abs(y).*besselk(nu,2/(pi*rh)*abs(y));
                % Compute the value at zero lag
                Cy(y==0)=s2;
            elseif nu==3/2
                % Second-order autoregressive
                Cy=s2*exp(-sqrt(6)/(pi*rh)*abs(y)).*(1+sqrt(6)/(pi*rh)*abs(y));
            elseif nu==5/2
                % Third-order autoregressive
                Cy=s2*exp(-sqrt(10)/(pi*rh)*abs(y)).*(1+sqrt(10)/(pi*rh)*...
                   abs(y)+10/(3*pi^2*rh^2)*abs(y).^2);            
            elseif mod(nu+0.5,1)==0
                % A general half-integer case beyond 1/2, 3/2, or 5/2, 
                % nu=n+1/2 for n = 1, 2, 3, ...
                k=0:nu-0.5;k=k(:); 
                n=nu-0.5;
                Cy=s2*exp(-2/(pi*rh)*abs(y)*sqrt(nu))*factorial(n)/factorial(2*n).*...
                   sum(factorial(n+k)./(factorial(k).*factorial(n-k)).*...
                   (4*sqrt(nu)*abs(y)/(pi*rh)).^(n-k),1); 
            else
                % However, if the nu provided to MATERNOSY is not one of the
                % special values of nu that we have considered so far, we should 
                % throw an error
                error('This is not a special case of nu. Ask for meth=1.')
            end
        end
    end
  % Or, calculate the derivatives?
  else 
      % Calculate the partial derivative wrt the requested parameter index
      % disp(sprintf('Calculating %i%s parameter derivative',dth,ith(dth)))
      % Abuse of Cy nomenclature, now you're getting a specific derivative
      if dth==1
          % The partial derivative of Cy with respect to the variance, dCyds2
          Cy=2^(1-nu)/gamma(nu)*argu.^nu.*besselk(nu,argu);
      elseif dth==2
          % The partial derivative of Cy with respect to the smoothness, 
          % dCydnu; the derivative of the gamma function is provided by Eq.
          % 6.3.1 of Abramowitz & Stegun (1965); 
          if nu==0
              % the partial derivative of the modified Bessel function of the
              % second kind with respect to order when the order is evaluated
              % at 0 is given by Eq. 9.6.46 of A&S, which we do not consider for
              % Matern
              dKdnuo=0;
          elseif mod(nu,1)==0
              % for integer order, the partial derivative of the modified
              % Bessel function of the second kind with respect to order is 
              % provided by Eq. 9.6.45 of A&S 
              syms k;
              dKdnuo=(factorial(nu).*(argu/2).^(-nu)/2).*...
                  vpasum((argu/2).^k.*besselk(k,argu)./...
                      ((nu-k).*factorial(k)),0,nu-1);
          elseif nu==1/2
              % for order of 1/2, we use Eq. 10.2.34 of A&S; note that
              % dCydnu(y=0,nu=1/2)=NaN because Ei(0)=-Inf; there also appear to
              % be numerical instability issues -- many NaNs and Infs for small
              % parameter values, e.g., th0=[1 0.5 1]; Olivia will study this
              % further. Default s2,rh for th0=[1e6 0.5 2e4] seems to behave. 
              warning('derivative wrt nu is still in development for nu given')
              dKdnuo=-sqrt(pi/2*argu).*ei(-2*argu).*exp(argu);
          %elseif mod(nu-1/2,1)==0
              % for half-integer order, Brychkov & Geddes (2005) presented a
              % solution for n-1/2 and n+1/2 order. We will use the n-1/2
              % order (Eq. 26) 
              %warning('derivative wrt nu is still in development for nu given')
              % I need to look into this expression further before we include
              % it.
              %n=nu+0.5;
              %argun=2*sqrt(n)/pi/rh*abs(y);
              %syms k p;
              %dKdnuo=(-1)^n*(pi/2)*(besseli(n-0.5,argun)+...
              %   besseli(0.5-n,argun)).*(coshint(2*argun)-sinhint(2*argun))+0.5*...
              %   vpasum((nchoosek(n,k)).*factorial(n-k-1).*(argun/2).^(k-n).*...
              %       besselk(k-0.5,argun),0,n-1)-...
              %   (-1)^n*sqrt(pi)*sum((-1).^k.*(nchoosek(n,k)).*2.^(k-1).*...
              %       (besseli(n-k-0.5,argun)+besseli(k-n+0.5,argun)).*...
              %        vpasum((nchoosek(k-1,p)).*factorial(k-p-1).*argun.^(p-k+0.5).*...
              %            besselk(p-0.5,2*argun),0,k-1),...
              %        1,n);
          else
              % for non-integer order, substitute DLMF Eq. 10.38.1 into 
              % Eq. 9.6.43 of A&S; we need to make an approximation to the 
              % infinite sum terms by choosing a large sumlim
              warning('derivative wrt nu is still in development for nu given')
              sumlim=1;
              syms k;
              dIdnup=besseli(nu,argu).*log(argu/2)-(argu/2).^nu.*...
                  vpasum((psi(nu+k+1)./gamma(nu+k+1).*...
                  (argu/2).^(k*2)./factorial(k)),0,sumlim);
              dIdnun=-besseli(-nu,argu).*log(argu/2)+(argu/2).^(-nu).*...
                  vpasum((psi(-nu+k+1)./gamma(-nu+k+1).*...
                  (argu/2).^(k*2)./factorial(k)),0,sumlim);
              dKdnuo=pi/2*csc(nu*pi)*(dIdnun-dIdnup)-...
                  pi*cot(nu*pi)*besselk(nu,argu);
              keyboard % solution blowing up for large y when th0=[1e6 1.5 1e6],
                       % NaNs and Inf when [1 1.5 1]?
              figure();plot(y(1:end-50),dKdnuo(1:end-50))
          end
          % Now put things together
          Cy=s2*2^(1-nu)/gamma(nu)*argu.^nu.*...
              (besselk(nu,argu).*(0.5+log(argu/2)-psi(nu))-...
              argu.^nu/(2*nu).*(argu.*besselk(nu-1,argu)+...
                                nu*besselk(nu,argu)).*double(dKdnuo));
      elseif dth==3
          % The partial derivative of Cy with respect to the range, dCydrho;
          % simplification of the derivative of the Bessel term with respect
          % to argument is made through Eq. 3.71.3 of Watson (1962)
          Cy=(s2/rh)*2^(1-nu)/gamma(nu)*argu.^(nu+1).*besselk(nu-1,argu);
      else
          error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
      end
    end
    % Serve output
    varns={Cy};
    varargout=varns(1:nargout);
elseif strcmp(y,'demo1')
    % The Fourier relation between the correlation and the covariance
    % is hard to verify... since we can never observe all lags on a
    % fine and large enough grid... and that is the point of estimating
    % parameters differently. See the comparison in SIMULOSL1('demo1')

    % Matern parameters... very sensitive to the last one
    th=[2000 0.5 2/3];
    % Grid size, also in physical units, keep it even for this example
    p.NyNx=[4300 5500]+randi(1000,[1 2])*(-1)^round(rand);
    p.NyNx=p.NyNx+mod(p.NyNx,2);
    lY=13; lX=24;
    p.dydx=[lY/(p.NyNx(1)-1) lX/(p.NyNx(2)-1)];
    % Wavenumber grid
    [k,kx,ky,dci]=knum2(p.NyNx,[lY lX]);
    % Space grid
    x=[-floor(p.NyNx(2)/2):1:+floor(p.NyNx(2)/2)-1]*p.dydx(2);
    y=[-floor(p.NyNx(1)/2):1:+floor(p.NyNx(1)/2)-1]*p.dydx(1);
    [X,Y]=meshgrid(x,y); yy=sqrt(X.^2+Y.^2);
    % Evaluate the Matern spectral covariance
    Sbb=v2s(maternos(k,th,[],2),p);
    % Evaluate the Matern correlation
    Cy=maternosy(yy,th,[],2); difer(Cy(dci(1),dci(2))-th(1))

    % Fourier transform the correlation to check its relation to the covariance
    Skk=fftshift(v2s(tospace(Cy,p),p));
    
    % Compare two profiles somehow
    m=median(Sbb(dci(1),:)./abs(Skk(dci(1),:)));
    m=max(max(Sbb))./max(max((Skk)));
    
    plot(log10(abs(Skk(dci(1),:))),'Color','r');
    hold on
    plot(log10(Sbb(dci(1),:))-log10(m),'LineWidth',2,'Color','m');
    hold off
    legend('FFT(MATERNOSY)','MATERNOS')
    axis tight
    grid on
    % Annotate 
    title(sprintf('p.NyNx = [%i %i]',p.NyNx))
    xlabel('wavenumber')
    ylabel('MATERNOS 2s FFT(MATERNOSY)')
elseif strcmp(y,'demo2')
    clf
    th0 = [1 1 10]; p = []; p.NyNx = [256 256]; p.dydx = [1 1]; 
    y = linspace(0,sqrt(prod(p.dydx))*sqrt(prod(p.NyNx)),1000);
    % Calculate Matern Covariance
    Cy = maternosy(y,th0,1);
    labs={'\sigma^2','nu','rho'};
    % Calculate Matern Covariance derivatives
    for index=1:3
        subplot(2,3,index)
        pc(index)=plot(y,Cy);
        ylim([0 th0(1)]+[-1 -1]*th0(1)/20)
        xlim([0 10*pi*th0(3)]+[-1 -1]*10*pi*th0(3)/20)
        title(sprintf('%s = %g','\nu',th0(2)))
        
        subplot(2,3,index+3)
        zdiff=maternosy(y,th0,index);
        pd(index)=plot(y,zdiff);
        keyboard
        ylim(halverange(zdiff,105,NaN))

        xlim([0 10*pi*th0(3)]+[-1 -1]*10*pi*th0(3)/20)
        title(labs{index})
    end
    keyboard
end
