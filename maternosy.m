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
% meth     1 straight up calculation (for finite nu) 
%          2 analytical simplification for special values of th(end-1)
%
% OUTPUT:
%
% CyordCydth    The spatial Matern covariance at all the requested lags, possibly
%               calculated using the Abramowitz & Stegun (1965) simplifications to
%               the Bessel function for half-integer values of nu, the infinite limit
%               of nu, OR returns the DERIVATIVE in the dth element of th, to feed 
%               into BLUROSY
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
% % C(0)=(2*pi)\int S(k)kdk
% p.NyNx=[bigN bigN];
% [k,dci,dcn,kx,ky]=knums(p); k=k(dci(1),dci(2):end);
% [sum(k(:).*maternos(k,th0,[],2))*(k(2)-k(1))*(2*pi) maternosy(0,th0) th0(1)]
% [sum(      maternos(k,th0,[],1))*(k(2)-k(1))*2      maternosy(0,th0) th0(1)]
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
% Calculate a partial derivative selected by the index provided as dth
%
% dths=[1 2 3]; dth=dths(randi(3)); dthlabs={'s2','nu','rh'};
% sprintf('Calculating the partial derivative w.r.t. %s',dthlabs{dth})
% dCydth=maternosy(y,th0,dth,1);
% plot(y,dCydth); % quick plot of partial derivative vs lag distance
%
% Demo1 provides a visual comparison of the general output of MATERNOS and
% MATERNOSY as Fourier pairs:
%
% maternosy('demo1')
%
% Demo2 is an earlier visualization of the partial derivatives (Keep?)
%
% maternosy('demo2')
%
% Checks the special cases
%
% maternosy('demo3')
%
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2025
% Last modified by olwalbert-at-princeton.edu, 04/26/2025

if ~isstr(y)
  % Defaults (avoiding DEFVAL to speed up)
  if ~exist('meth','var') | isempty(meth); meth=1; end
  if ~exist('dth','var'); dth=[]; end

  % The Matern parameters are always the last three elements of the TH input 
  s2=th(end-2);
  nu=th(end-1);
  rh=th(end  );
  % The argument, make sure it is a distance
  argu=2*sqrt(nu)/pi/rh*abs(y);

  % Check whether we are asking for the case that nu approaches
  % infinity; if so, by the straight-up method Cy would have Inf/NaN entries
  if isinf(nu)
    % Squared exponential 
    % This analytic form of the isotropic Matern covariance 
    % was solved from the inverse Fourier transform of the limit as
    % nu approaches infinity of the spectral density (see MATERNOS
    % for details) with use of Eq. 3.323.2 of Gradshteyn & Ryzhik (1980).
    if isempty(dth)
        % Calculate the spatial covariance
        Cy=s2/(pi*rh)*exp(-abs(y).^2/(pi^2*rh^2));
    elseif dth==1
        % Calculate the partial derivative with respect to s2 
        Cy=exp(-abs(y).^2/(pi^2*rh^2));
    elseif dth==2
        % Calculate the partial derivative with respect to nu
        Cy=zeros(size(y));
    elseif dth==3
        % Calculate the partial derivative with respect to rho
        Cy=s2/rh*exp(-y.^2/(pi^2*rh^2))*(1/(pi*rh)).*(2*y.^2/(pi^2*rh^2)-1);
    else
        error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
    end
  else
    % Switch the calculation method
    switch meth
      case 1
        if isempty(dth)
            % Calculate the spatial covariance
            % The evaluation
            Cy=2^(1-nu)*s2/gamma(nu)*argu.^nu.*besselk(nu,argu);
            % Supply the smallest arguments
            Cy(y==0)=s2;
        else
          % Calculate the partial derivative wrt the requested parameter index
          % disp(sprintf('Calculating %i%s parameter derivative',dth,ith(dth)))
          % Abuse of Cy nomenclature, now you're getting a specific derivative
          if dth==1
              % The partial derivative of Cy with respect to the variance, dCyds2
              Cy=2^(1-nu)/gamma(nu)*argu.^nu.*besselk(nu,argu);
              % Supply the smallest arguments from the partial derivative for the 
              % small argument limit of the Bessel function, see Eq. 9.6.9 of A&S
              Cy(y==0)=1;
          elseif dth==2
            % The partial derivative of Cy with respect to the smoothness, 
            % dCydnu; the derivative of the gamma function is provided by Eq.
            % 6.3.1 of Abramowitz & Stegun (1965); 

            % We are still experiencing numerical issues. Turn off the
            % analytical solution for now and calculate via finite difference.
            anlytcl=0;
            if anlytcl
              if nu==0
                  % The partial derivative of the modified Bessel function of the
                  % second kind with respect to order evaluated for order==0
                  % is given by Eq. 9.6.46 of A&S; (we do not consider this case) 
                  dKdnuo=0;
              elseif (mod(nu,1)~=1)&((ceil(nu)-nu)>=0.05)
                  % We calculate the derivative of the Bessel function term 
                  % following in the steps of Brychkov (2016) using Mathematica,
                  % with the added complication of nu appearing in the Matern 
                  % argument; this derivative is valid for positive, non-integer
                  % values of nu. We are trying alternatives to MATLAB's
                  % HYPERGEOM, which relies on conversion to symbolic values
                  % (very, very slow)
                  
                  % Pre-compute common terms to save time
                  BKNZ=besselk(nu,argu);
                  % BINPZ may be the source of the numerical problems for large
                  % ARGU
                  BINPZ=besseli(nu,argu); % positive order
                  CSCNUPI=csc(nu*pi);
                  ARGUSQ=argu.^2;
                  HFARGUSQ=(argu/2).^2;

                  % Most recent attempt using genHyper from file exchange;
                  % GENHYPER does not take vector arguments
                  NumWorkers=feature('numcores')-1;
                  if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end
                  parfor ind=1:numel(argu)
                   GENHYPER1(ind)=genHyper([nu,nu+1/2],[nu+1,nu+1,2*nu+1],ARGUSQ(ind));
                   GENHYPER2(ind)=genHyper([1/2-nu,-nu],[1-2*nu,1-nu,1-nu],ARGUSQ(ind));
                   GENHYPER3(ind)=genHyper([1,1,3/2],[2,2,2-nu,2+nu],ARGUSQ(ind));
                  end
                  t1=(pi/2.*CSCNUPI).^2*(HFARGUSQ.^(nu)./(gamma(1+nu).^2).*...
                     besseli(-nu,argu).*GENHYPER1-HFARGUSQ.^(-nu)/(gamma(1-nu).^2).*BINPZ.*...
                     GENHYPER2);
                  t2=1/(2*nu).*(argu.*besselk(nu-1,argu)+(nu-1).*BKNZ-...
                     (pi+pi^2.*nu.*cot(pi*nu)).*CSCNUPI.*BINPZ);
                  t3=(BKNZ+pi*CSCNUPI.*BINPZ).*(HFARGUSQ./(nu.^2-1).*...
                     GENHYPER3-log(argu)+psi(nu));
                  % Original, using MATLAB's hypergeom
                  % t1=(pi/2.*CSCNUPI).^2*(HFARGUSQ.^(nu)./(gamma(1+nu).^2).*...
                  %    besseli(-nu,argu).*hypergeom([nu,nu+1/2],[nu+1,nu+1,2*nu+1],ARGUSQ)-...
                  %    HFARGUSQ.^(-nu)/(gamma(1-nu).^2).*BINPZ.*...
                  %    hypergeom([1/2-nu,-nu],[1-2*nu,1-nu,1-nu],ARGUSQ));
                  % t2=1/(2*nu).*(argu.*besselk(nu-1,argu)+(nu-1).*BKNZ-...
                  %    (pi+pi^2.*nu.*cot(pi*nu)).*CSCNUPI.*BINPZ);
                  % t3=(BKNZ+pi*CSCNUPI.*BINPZ).*(HFARGUSQ./(nu.^2-1).*...
                  %    hypergeom([1,1,3/2],[2,2,2-nu,2+nu],ARGUSQ)-log(argu)+psi(nu));
                  dKdnuo=t1-t2+t3;
              else %mod(nu,1)==0
                  % For integer order, we must calculate the limit from the
                  % non-integer order form of the modified Bessel function of the 
                  % second kind; note that the singularity that arises at integer
                  % order is due to the cosecant term in the derivative

                  % For now, let's do this through averaging on either side of
                  % the singularity (is there a nicer way?)
                  pert=0.01;
                  thleft=th;thleft(2)=nu-pert;
                  dCydnuleft=maternosy(y,thleft,2);
                  thright=th;thright(2)=nu+pert;  
                  dCydnuright=maternosy(y,thright,2);

                  dKdnuo=0.5*(dKdnuoleft+dKdnuoright);
              end
              % Now put things together
              Cy=s2*2^(1-nu)/gamma(nu)*argu.^nu.*...
                  (besselk(nu,argu).*(0.5+log(argu/2)-psi(nu))+...
                  double(dKdnuo));
              % Supply the smallest arguments from the partial derivative for the 
              % small argument limit of the Bessel function, see Eq. 9.6.9 of A&S
              Cy(y==0)=0;
            else
              % Central difference approach for calculating the partial
              % derivative that we default to.
              % Set the step-size; this first value works well in practice
              h=eps*5e5/2;
              % Press+2007 5.7.8 list the optimal step-size
              h=eps^(1/3)*th(2); 
              % We know there are discontinuities at integer values, so if th(2)
              % +/- h is close to an integer, double h
              while abs(mod(th(2)+h,1))<1e-10 | abs(mod(th(2)-h,1))<1e-10
                 h=2*h;
              end
              % Parameters evaluated for step-size around th
              thl=th;thl(2)=thl(2)-h;
              thr=th;thr(2)=thr(2)+h;
              % Autocovariance evaluated around th(2)
              Cy=(maternosy(y,thr)-maternosy(y,thl))/(2*h);
              % Press+2007 eq 5.7.9: error from roundoff and truncation is 
              % approx eps^(2/3) for this finite difference formulation and
              % step-size

             % Another alternative is chebyshev polynomial approximation for Cy
             % as a function of nu that can be differentiated
             % k=20; t=augknt(minmax(y),k); x=chbpnt(t,k);
             % %sp=spapi(t,x,maternosy(x,th0,2,1));
             % % plot(y,fnval(sp,y)); hold on; plot(y,maternosy(y,th0,2,1),'--')
             % sp=spapi(t,x,maternosy(x,th0,[],1));
             % Cy=fnval(sp,y);
             % % evaluate:
             % plot(y,Cy); hold on; plot(y,maternosy(y,th0,[],1),'--')
            end
          elseif dth==3
              % The partial derivative of Cy with respect to the range, dCydrho;
              % simplification of the derivative of the Bessel term with respect
              % to argument is made through Eq. 3.71.3 of Watson (1962)
              Cy=(s2/rh)*2^(1-nu)/gamma(nu)*argu.^(nu+1).*besselk(nu-1,argu);
              % Supply the smallest arguments from the partial derivative for the 
              % small argument limit of the Bessel function, see Eq. 9.6.9 of A&S
              Cy(y==0)=0;
          else
              error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
          end
        end
      case 2
      % By selecting the second method of calculation, we are seeking to 
      % evaluate Cy from the simplified analytic expression of the 
      % isotropic Matern covariance for a special value of nu. The 
      % following analytic forms of half-integer nu are calculated from 
      % substitution of Eqs 10.1.9 and 10.2.15 of Abramowitz & Stegun 
      % (1965) for modified Bessel functions of half-integer orders.
      % See also 10.1093/biomet/93.4.989 for comparison.
        if nu==1/3
          if isempty(dth)
            % Calculate the spatial covariance
            % von Karman 
            Cy=s2*2^(2/3)/gamma(nu).*(2*sqrt(3)/(3*pi*rh)*y).^...
               (nu).*besselk(nu,2*sqrt(3)/(3*pi*rh)*y);
            % Compute the value at zero lag
            Cy(y==0)=s2;
          elseif dth==1
            %Calculate the partial derivative with respect to s2
            Cy=2^(2/3)/gamma(nu).*(2*sqrt(3)/(3*pi*rh)*y).^...
               (nu).*besselk(nu,2*sqrt(3)/(3*pi*rh)*y);
            % Compute the value at zero lag
            Cy(y==0)=1;
          elseif dth==2
            %Calculate the partial derivative with respect to nu
            Cy=zeros(size(y));
          elseif dth==3
            %Calculate the partial derivative with respect to rho
            Cy=s2/rh*2^(2/3)/gamma(1/3).*(2*sqrt(3)/(3*pi*rh)*y).^(4/3).*...
               besselk(2/3,2*sqrt(3)/(3*pi*rh)*y);
            % Compute the value at zero lag
            Cy(y==0)=s2;
          else
            error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
          end
        elseif nu==1/2
          if isempty(dth)
            % Calculate the spatial covariance
            % Exponential
            Cy=s2*exp(-sqrt(2)/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=s2;
          elseif dth==1
            %Calculate the partial derivative with respect to s2
            Cy=exp(-sqrt(2)/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=1;
          elseif dth==2
            %Calculate the partial derivative with respect to nu
            Cy=zeros(size(y));
          elseif dth==3
            %Calculate the partial derivative with respect to rho
            Cy=s2/rh*exp(-sqrt(2)/(pi*rh)*abs(y)).*(sqrt(2).*y/(pi*rh));
            % Compute the value at zero lag
            %%% OLW TODO: Confirm that the following is correct
            % Cy(y==0)=s2;
          else
            error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
          end
        elseif nu==1
          if isempty(dth)
            % Calculate the spatial covariance
            % Whittle
            Cy=s2*2/(pi*rh)*abs(y).*besselk(nu,2/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=s2;
          elseif dth==1
            %Calculate the partial derivative with respect to s2
            Cy=2/(pi*rh)*abs(y).*besselk(nu,2/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=1;
          elseif dth==2
            %Calculate the partial derivative with respect to nu
            Cy=zeros(size(y));
          elseif dth==3
            %Calculate the partial derivative with respect to rho
            Cy=s2/rh*(2/(pi*rh)*abs(y)).^2.*besselk(0,2/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=s2;
          else
            error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
          end
        elseif nu==3/2
          if isempty(dth)
            % Calculate the spatial covariance
            % Second-order autoregressive
            Cy=s2*exp(-sqrt(6)/(pi*rh)*abs(y)).*(1+sqrt(6)/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=s2;
          elseif dth==1
            %Calculate the partial derivative with respect to s2
            Cy=exp(-sqrt(6)/(pi*rh)*abs(y)).*(1+sqrt(6)/(pi*rh)*abs(y));
            % Compute the value at zero lag
            Cy(y==0)=1;
          elseif dth==2
            %Calculate the partial derivative with respect to nu
            Cy=zeros(size(y));
          elseif dth==3
            %Calculate the partial derivative with respect to rho
            Cy=s2/rh*exp(-sqrt(6)/(pi*rh)*abs(y)).*(sqrt(6)/(pi*rh)*abs(y)).^2;
            % Compute the value at zero lag
            Cy(y==0)=s2;
          else
            error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
          end
        elseif nu==5/2
          if isempty(dth)
            % Calculate the spatial covariance
            % Third-order autoregressive
            Cy=s2*exp(-sqrt(10)/(pi*rh)*y).*(1+sqrt(10)/(pi*rh)*y+...
               10/(3*pi^2*rh^2)*y.^2);            
            % Compute the value at zero lag
            Cy(y==0)=s2;
          elseif dth==1
            %Calculate the partial derivative with respect to s2
            Cy=exp(-sqrt(10)/(pi*rh)*y).*(1+sqrt(10)/(pi*rh)*y+...
               10/(3*pi^2*rh^2)*y.^2);            
            % Compute the value at zero lag
            Cy(y==0)=1;
          elseif dth==2
            %Calculate the partial derivative with respect to nu
            Cy=zeros(size(y));
          elseif dth==3
            %Calculate the partial derivative with respect to rho
            Cy=s2/3/rh*exp(-sqrt(10)/(pi*rh)*y).*(sqrt(10)/(pi*rh)*y).^2.*...
                (1+sqrt(10)/(pi*rh)*y);
            % Compute the value at zero lag
            Cy(y==0)=s2;
          else
            error('Not a valid partial derivative index. Ask for dth=1,2,or 3.')
          end
        elseif mod(nu+0.5,1)==0
          if isempty(dth)
            % Calculate the spatial covariance
            % A general half-integer case beyond 1/2, 3/2, or 5/2, 
            % nu=n+1/2 for n = 1, 2, 3, ...
            ydim=size(y);                                                       
            % If y is a row-vector, make it a column vector                     
            if ydim(1)==1 & ydim(2)>1; y=y(:); end                              
            % If y is 2+-dimensional, push k to the next dimension        
            k=0:nu-0.5; k=reshape(k,[repelem(1,numel(size(y)>1)) numel(k)]);
            n=nu-0.5;
            Cy=s2*exp(-2/(pi*rh)*abs(y)*sqrt(nu))*factorial(n)/factorial(2*n).*...
               sum(factorial(n+k)./(factorial(k).*factorial(n-k)).*...          
                   (4*sqrt(nu)*abs(y)/(pi*rh)).^(n-k),numel(size(y)>1)+1);          
            Cy=reshape(Cy,ydim);           
            % Compute the value at zero lag
            Cy(y==0)=s2;
          else
            warning(append('Partial derivatives for the half-integer special',...
              'case are in progress. Returning the requested partial',...
              'derivative for the general Matern form instead'))
            Cy=maternosy(y,th,dth,1);
          end
        else
          % However, if the nu provided to MATERNOSY is not one of the
          % special values of nu that we have considered so far, we should 
          % throw an error
          warning(append('This is not a special case of nu. Returning the',...
                         'requested autocovariance for the general Matern form instead'))
          Cy=maternosy(y,th,dth,1);
        end
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

    % Physical length
    lY=13; lX=24;

    % Matern parameters... think about last one in terms of sqrt(lX^2+lY^2)
    th=[2000 4.5 2/6];

    % Grid size, also in physical units, keep it even for this example
    p.NyNx=[4300 5500]; %+randi(1000,[1 2])*(-1)^round(rand);
    p.NyNx=p.NyNx+mod(p.NyNx,2);
    p.dydx=[lY/(p.NyNx(1)-1) lX/(p.NyNx(2)-1)];
    % Wavenumber grid
    [k,kx,ky,dci]=knum2(p.NyNx,[lY lX]);
    % Space grid
    x=[-floor(p.NyNx(2)/2):1:+floor(p.NyNx(2)/2)-1]*p.dydx(2);
    y=[-floor(p.NyNx(1)/2):1:+floor(p.NyNx(1)/2)-1]*p.dydx(1);
    [X,Y]=meshgrid(x,y);
    yy=sqrt(X.^2+Y.^2);
    keyboard
    % Evaluate the Matern spectral covariance
    methS=randi(2); methS=2;
    Sbb=v2s(maternos(k,th,[],2,methS),p);
    % Evaluate the Matern correlation
    methC=randi(2); methC=2;
    Cy=maternosy(yy,th,[],methC);
    % Check the variance
    difer(Cy(dci(1),dci(2))-th(1))

    % ONE-DIMENSIONAL 
    subplot(211)
    % Do the one-dimensional profile through the 2-D thing
    % bb1=Sbb(dci(1),dci(2):end);
    % This is not the same as the 1-dimensional Matern
    Sbb1=maternos(k(dci(1),dci(2):end),th,[],1);
    % This is the one-dimensional correlation
    Cy1=Cy(dci(1),:);
    % Ignore the hopefully tiny imaginary parts
    %fCy1=realize(fft(Cy1));
    ffCy1=fftshift(realize(fft(fftshift(Cy1))));
    Skk1=ffCy1(dci(2):end)*p.dydx(2)/(2*pi);

    semilogx(Sbb1,'o-')
    hold on
    semilogx(Skk1,'+')
    hold off
    grid on
    legend('FFT(MATERNOSY)','MATERNOS 1D')
    axis tight

    % TWO-DIMENSIONAL
    % Fourier transform the correlation to check its relation to the covariance
    % Cy is already in space; the following two are identical
    %Skk=v2s(tospec(Cy,p),p)*(sqrt(prod(p.dydx))*sqrt(prod(p.NyNx)))/(2*pi)^2;
    Skk=fftshift(fft2(Cy)*prod(p.dydx)/(2*pi)^2);

    subplot(212)
    %plot(log10(Sbb(dci(1),:))-log10(m),'LineWidth',2,'Color','m');
    %plot(log10(kx(dci(2):end)),(Sbb(dci(1),dci(2):end)),'LineWidth',2,'Color','m');
    semilogx(Sbb(dci(1),dci(2):end),'o-');
    hold on
    %plot(log10(abs(Skk(dci(1),:))),'Color','r');
    %plot(log10(kx(dci(2):end)),(abs(Skk(dci(1),dci(2):end))),'Color','r');
    semilogx(abs(Skk(dci(1),dci(2):end)),'+');
    hold off
    grid on
    legend('FFT(MATERNOSY)','MATERNOS 2D')
    axis tight

    % Annotate 
    title(sprintf('p.NyNx = [%i %i]',p.NyNx))
    xlabel('wavenumber')
    ylabel('MATERNOS 2s FFT(MATERNOSY)')
elseif strcmp(y,'demo2')
    clf
    th0 = [1 4 20]; p = []; p.NyNx = [256 256]; p.dydx = [1 1];

    y = linspace(0,sqrt(prod(p.dydx)*prod(p.NyNx)),1000);
    xels=[0 min(6*pi*th0(3),max(y))];
    xels=xels+[-1 1]*range(xels)/20;
    
    % Calculate Matern Covariance
    labs={'\sigma^2','\nu','\rho'};
    % Calculate Matern Covariance derivatives
    for index=1:3
        ah1(index)=subplot(2,3,index);
        pc(index)=plot(y,maternosy(y,th0,1),'k','LineWidth',1); hold on
        % Let's add a few neighbors

        % Percentage perturbation in each of the three parameters 
        pc1(index)=plot(y,maternosy(y,th0.*(1+rindeks(eye(3),index)/20)),'r');
        pc2(index)=plot(y,maternosy(y,th0.*(1-rindeks(eye(3),index)/20)),'b'); hold off

        ylim([0 th0(1)]+[-1 1]*th0(1)/20)
        %        ylim([-0.1 0.1])
        xlim(xels)
        if index==2
            title(sprintf('%s = %g  | %s = %g | %s = %g',...
                          '\sigma^2',th0(1),'\nu',th0(2),'\rho',th0(3)))
        end
        grid on
        ylabel('covariance C(y)')
        xlabel('lag y')
                
        ah2(index)=subplot(2,3,index+3);
        % Calculate derivative
        zdiff=maternosy(y,th0,index);
        pd(index)=plot(y,zdiff,'k');
        ylim(halverange(zdiff,105,NaN))

        % if index==2
        %     hold on
        %     plot([1 1]*pi*th0(3)/2/sqrt(th0(2)),ylim)
        %     hold off
        % end
        grid on
        ylabel(sprintf('covariance derivative dC(y)/d%s',labs{index}))
        xlabel('lag y')
        xlim(xels)
    end
    longticks([ah1 ah2],2)
elseif strcmp(y,'demo3')
    s2=rand*100;
    nus=[1/3 1/2 1 3/2 5/2 3.5 4.5];
    rho=rand*100;
    for index=1:length(nus)
        th=[s2 nus(index) rho];
        % Just create some lags out to some number of times the correlation length 
        y=linspace(0,pi*th(3)*4,100);
        Cy1=maternosy(y,th,[],1);
        Cy2=maternosy(y,th,[],2);
        difer(Cy1-Cy2)
        % Also test MATERNOS
        % Just create some wavenumbers that are appropriate
        k=sort([2*pi./y(end:-1:2) 0]);
        Sk1=maternos(k,th,[],2,1);
        Sk2=maternos(k,th,[],2,2);
        difer(Sk1-Sk2,9)
    end
end
