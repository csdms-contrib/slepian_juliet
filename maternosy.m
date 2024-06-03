function varargout=maternosy(y,th,meth,dth)
% varargout=MATERNOSY(y,th,meth,dth)
%
% Calculates d-dimensional isotropic Matern correlation, which is
% independent of d, also the counterpart to the spectral covariance.
% See Olhede & Simons (2013), doi: 10.1093/gji/ggt056.x, eq. (72)
%
% INPUT:
%
% y        Lag parameter, the distance between spatial positions,
%          e.g. from XXPDIST or SSPDIST [m], can be any dimension
% th       The unscaled parameter vector, with, in the last three slots: 
%          s2    The first Matern parameter [variance in units^2]
%          nu    The second Matern parameter [differentiability]
%          rh    The third Matern parameter [range in units]
% meth     1 straight up calculation
%          2 analytical simplification for special values of th(end-1)
% dth      1, 2, or 3 specifies which element of th gets differentiated
%
% OUTPUT:
%
% Cy       The spatial Matern covariance at all the requested lags, possibly
%          calculated using the Abramowitz & Stegun (1965) simplifications to
%          the Bessel function for half-integer values of nu, OR
%          Returns the DERIVATIVE in the dth element of th, to feed into BLUROSY
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
% nu=[1/3 1/2 1 3/2 5/2]; th0(2)=nu(randi(length(nu))); th0(2)
% Cy=maternosy(y,th0,1); Cy2=maternosy(y,th0,2); difer(Cy/(Cy(1))-Cy2/Cy2(1))
%
% th0(2)=Inf; Cy=maternosy(y,th0);
%
% maternosy('demo1')
%
% Last modified by fjsimons-at-alum.mit.edu, 06/03/2024
% Last modified by olwalbert-at-princeton.edu, 06/03/2024

if ~isstr(y)
    % Defaults
    defval('meth',1)
    defval('dth',[])
    % These are always the last three elements of the input 
    s2=th(end-2);
    nu=th(end-1);
    rh=th(end  );
    % The argument, make sure it is a distance
    argu=2*sqrt(nu)/pi/rh*abs(y);
    % Check whether we are asking for the case that nu approaches
    % infinity; if so, by the regular method Cy would have Inf/NaN entries
    if isinf(nu)
        % Squared exponential 
        % This analytic form of the isotropic Matern covariance 
        % was solved from the inverse Fourier transform of the limit as
        % nu approaches infinity of the spectral density (see MATERNOS
        % for details). Eq. 3.323.2 of Gradshteyn & Ryzhik (1980) was
        % applied. (Consistency with BLUROSY demo?)
        Cy=s2/(pi*rh)*exp(-abs(y).^2/(pi^2*rh^2));
    else
        % Switch the calculation method
        switch meth
          case 1
            % The evaluation, noting that nu=Inf will be overwritten properly
            Cy=2^(1-nu)*s2/gamma(nu)*argu.^nu.*besselk(nu,argu);
            % Supply the smallest arguments
            Cy(y==0)=s2;
          case 2
            % If the number of outputs requested is greater than 1, we
            % are seeking to evaluate Cy2 from the simplified analytic
            % expression of the isotropic Matern covariance for a special 
            % value of nu. The following analytic forms of half-integer
            % nu are calculated from substitution of Eqs 10.1.9 and 10.2.15
            % of Abramowitz & Stegun (1965) for modified
            % Bessel functions of half-integer orders.
            % See also 10.1093/biomet/93.4.989
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
            else
                % However, if the nu provided to MATERNOSY is not one of the
                % five special values of nu, we should throw an error
                error('This is not a special case of nu. Ask for meth=1.')
            end
        end
    end
    % Calculate the derivatives?
    if ~isempty(dth)
        % Calculate derivatives in the requested coordinate
        % disp(sprintf('Calculating %i%s parameter derivative',dth,ith(dth)))
        % Abuse of nomenclature, now you're getting a specific derivative
        Cy='beautiful';
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
    Sbb=v2s(maternos(k,th,2),p);
    % Evaluate the Matern correlation
    Cy=maternosy(yy,th,2); difer(Cy(dci(1),dci(2))-th(1))

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
    ylabel('MATERNOS vs FFT(MATERNOSY)')
end
