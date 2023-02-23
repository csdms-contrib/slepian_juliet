function varargout=maternosy(y,th,varargin)
% varargout=MATERNOSY(y,th,d)
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
% d        The dimensionality. Not used, for symmetry with MATERNOS only.
%
% OUTPUT:
%
% Cy       The spatial Matern covariance at all the requested lags
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
% [sum(y.*maternosy(y,th0))*(y(2)-y(1))/(2*pi) maternos(0,th0)]
%
% maternosy('demo1')
% 
% Last modified by fjsimons-at-alum.mit.edu, 02/23/2023

if ~isstr(y)
    % These are always the last three elements of the input 
    s2=th(end-2);
    nu=th(end-1);
    rh=th(end  );
    
    % The argument, make sure it is a distance
    argu=2*sqrt(nu)/pi/rh*abs(y);
    % The evaluation
    Cy=2^(1-nu)*s2/gamma(nu)*argu.^nu.*besselk(nu,argu);
    % Supply the smallest arguments
    Cy(y==0)=s2;
    % Optional output
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
    axis tight
    grid on
    % Annotate 
    title(sprintf('p.NyNx = [%i %i]',p.NyNx))
    xlabel('wavenumber')
    ylabel('MATERNOS vs FFT(MATERNOSY)')
end
