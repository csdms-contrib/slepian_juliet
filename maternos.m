function varargout=maternos(k,th,dth,d,meth)
% [SkordSkdth,k]=MATERNOS(k,th,dth,d,meth)
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
% dth      When empty, the spectral density is calculated; otherwise,
%          1, 2, or 3 specifies which element of th gets differentiated
% d        The dimensionality [default is 2]
% meth     1 straight up calculation
%          2 analytical simplification for special values of th(end-1)
%
% OUTPUT:
%
% SkordSkdth  A column vector with all the wavenumbers unwrapped, possibly
%             calculated using a simplification of nu, OR returns the DERIVATIVE in
%             the dth element of th; simplifications and derivatives are included 
%             to mirror MATERNOSY, however, these do not feed into another program
% k           Return the wavenumber grid, which was an input, so it works like BLUROS
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
% imagesc(log10(v2s(Sk,p))); axis image
%
% Confirmation that the special case simplifications are equivalent to the
% straight calculation:
%
% nus=[1/3,1,1/2,3/2,5/2,7/2,15/2]; th=[1 nus(randi(length(nus),1)) 1];
% Sk=maternos(k,th,[],2,1); Sk2=maternos(k,th,[],2,2); diferm(Sk,Sk2)
%
% Calculate and visualize the partial derivatives: 
%
% dth=randi(3,1); dSkdth=maternos(k,th,dth,2,1); 
% figure(); imagesc(log10(v2s(dSkdth,p))); axis image
%
% Demo1 provides a visual comparison of the spectral density, its partial
% derivatives with respect to the Matern parameters, and their central
% difference approximations; the scaling for s2 and rho will be improved
%
% maternos('demo1',[],dth)
%
% Last modified by fjsimons-at-alum.mit.edu, 08/25/2026
% Last modified by olwalbert-at-princeton.edu, 11/02/2026

if ~isstr(k)
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
                      Sk=s2*2^(2/3)*pi*rh^2./(4+3*(pi*rh*k(:)).^2).^(4/3); 
                    else
                      Sk=gamma(1/3+d/2)/gamma(1/3)*s2*...
                          (4./(4+3*(pi*rh*k(:)).^2)).^(1/3).*...
                          (3*pi*rh^2./(4+3*(pi*rh*k(:)).^2)).^(d/2);
                    end
                elseif nu==1/2
                    % Exponential
                    if d==2
                        Sk=s2*pi*rh^2./(sqrt(2)*(2+(pi*rh*k(:)).^2).^(3/2));
                    else
                        Sk=gamma((d+1)/2)*s2*(2/(pi^(1-d)))^(1/2)*rh^d*...
                            (2+(pi*rh*k(:)).^2).^(-(d+1)/2);
                    end
                elseif nu==1
                    % Whittle
                    if d==2
                      Sk=s2*4*pi*rh^2./(4+(pi*rh*k(:)).^2).^2; 
                    else
                      Sk=gamma(1+d/2)*s2/pi^(d/2)*(2/(pi*rh))^2*...
                          ((pi*rh)^2./(4+(pi*rh*k(:)).^2)).^(1+d/2);
                    end
                elseif nu==3/2
                    % Second-order autoregressive
                    if d==2
                      Sk=9*sqrt(6)*s2*pi*rh^2./(6+(pi*rh*k(:)).^2).^(5/2);
                    else
                      Sk=gamma((d+3)/2)*12*(6*pi)^(1/2)*s2*rh^2./...
                          (6+(pi*rh*k(:)).^2).^(5/2);
                    end
                elseif nu==5/2
                    % Third-order autoregressive
                    if d==2
                      Sk=250*sqrt(10)*s2*pi*rh^2./(10+(pi*rh*k(:)).^2).^(7/2);
                    else
                      Sk=gamma((d+5)/2)*s2*400/3*(10)^(1/2)*pi^(d/2-1/2)*rh^d*...
                          (10+(pi*rh*k(:)).^2).^(-(d+5)/2);
                    end
                elseif mod(nu-1/2,1)==0
                    % A general half-integer case beyond 1/2, 3/2, or 5/2,
                    % nu=n+1/2 for n = 1, 2, 3, ...
                    n=nu-1/2;
                    if d==2
                      Sk=s2*pi*rh^2/4*((4*n+2)./(4*n+2+(pi*rh*k(:)).^2)).^(n+3/2);
                    else
                      Sk=gamma(n+(d+1)/2)/gamma(n+1/2)*s2*((4*n+2)./...
                          (4*n+2+(pi*rh*k(:)).^2)).^(n+1/2).*...
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
     avark=4*nu/pi^2/rh^2+k(:).^2;
     if dth==1
         % The partial derivative of Sk with respect to the variance, dSkds2
         Sk=pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2);
     elseif dth==2
         % The partial derivative of Sk with respect to the smoothness, dSkdnu; 
         % see Eq. 6.3.1 of Abramowitz & Stegun (1965) for the derivative of the 
         % Gamma function and 
         Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2).*...
            (psi(nu+d/2)-psi(nu)+log(nu*4/(pi*rh)^(2))-log(avark)-...
            (4*nu+2*d)./(4*nu+(pi*rh*k(:)).^2)+1);
     elseif dth==3
         % The partial derivative of Sk with respect to the range, dSkdrh
         Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rh)^(2*nu).*avark.^(-nu-d/2-1)*...
            (2*nu/rh).*(2*d/(pi*rh)^2-k(:).^2);
     else
         error('Not a valid partial derivative index. Ask for dth=1,2,or 3.') 
     end
  end
  % Serve output
  varns={Sk};
  varargout=varns(1:nargout);
elseif strcmp(k,'demo1')
  % Creates a plot showing the isotropic, stationary Matern spectrum, 
  % its partial derivative with respect to a Matern parameter selected by the 
  % dth index in context of the wave-vector and the parameter 

  % The Matern parameters
  defval('th',[1e6 1.5 2e4])
  % Matern parameter index of the parameter that will be allowed to vary
  defval('dth',2)
  % The grid
  defval('p',[])
  if isempty(p)                                                                   
    p.NyNx=[128 128];                                                           
    p.dydx=[1e3 1e3];                                                           
    p.blurs=Inf;                                                                
  end  
  % How many dimensions?
  dim=2; 
  % How many curves to plot?
  numcs=3;
  % Set the x-, and when relevant, the y-axis to log-scale?
  ll=1;

  % For wavevector dependent subfigure, define a dense wavevector and devise
  % indices for numcs values of th(dth) that are of interest; for parameter
  % dependent subfigure, define a dense parameter vector for th(dth) and devise
  % indices for wavenumbers of interest
  
  % Wavenumber axes, 2D axes flattened to 1D radial, inspired from EGGERS2
  mfromkm=1000;
  Nx=p.NyNx(2);
  dx=p.dydx(2);
  % Wavenumber limits [rad/m]
  [k,~,~,dci]=knum2([1 Nx],[1 (Nx-1)*dx]);
  k=k(dci(2):end);
  % Finer grid for S(k)
  k=unique([linspace(k(1),k(2),50)...
      linspace(k(2),k(3),50)...
      linspace(k(3),k(4),50)...
      k(4:end)]);
  % Spatial axes, for distance [m]
  x=linspace(0,(Nx-1)*dx,100);
  % New range and ticks for Ky
  newr=[-0.1 1.1];
  newy=[0 1/3 2/3 1];
  % Major tick wavelengths [km]
  xtil=[100 1000 10000];
  % Wavenumber limits [rad/km]
  klims=[k(2) max(k)]*mfromkm;
  
  % Labels in order of dth indices
  labs={'\sigma^2','\nu','\rho'};
  axlabs={'variance, \sigma^2','smoothness, \nu','range, \rho'};
  
  % Matern parameters of interest
  thalld={[th(1)/2:th(1)/50:th(1)*3],[0.05:0.1:20],[th(3)/2:th(3)/50:th(3)*3]};
  thsub={[2*th(1)/3 th(1) 2*th(1)],[1/2 3/2 5/2],[2*th(3)/3 th(3) th(3)*2]};
  thall{dth}=[thalld{dth} thsub{dth}];
  % The fluctuation scale [m^4], following from EGGERS2; will need to be          
  % recalculated if dth==1 or dth==3                                              
  switch dth                                                                      
      case 1                                                                      
          S0d=thalld{dth}.*pi.*th(3).^2/4;                                        
          S0=thall{dth}.*pi.*th(3).^2/4;                                          
      case 2                                                                      
          S0d=th(1).*pi.*th(3).^2/4;                                              
          S0=S0d;                                                                 
      case 3                                                                      
          S0d=th(1).*pi.*thalld{dth}.^2/4;                                        
          S0=th(1).*pi.*thall{dth}.^2/4;                                          
  end                                                                             
  [thall{dth},temp]=sort(thall{dth});[~,thind]=sort(temp);                        
  thind=thind(end-numcs+1:end);                                                    
  N0=length(thalld{dth});                                                         
  N1=length(thsub{dth});                                                          
  N2=length(thall{dth});                                                          
                                                                                  
  % Wavenumbers of interest                                                       
  knd=[floor(length(k)/20); floor(length(k)/9); floor(length(k)/5)];              
  ksub=k(knd)';               
  
  % Preallocate space for the periodograms, their derivatives, and their finite
  % difference approximations
  Sks=zeros(N0,length(k));
  Sk1s=zeros(N1,length(k));
  Sk2s=zeros(N0,length(ksub));
  dSk1s=zeros(N1,length(k));
  dSk2s=zeros(N0,length(ksub));
  cdiffSks=zeros(length(k),N2);
  % Preallocate space for the L2 norms of ratios between periodograms of
  % increasing th
  rSk=zeros(N2-1,1);
  rdSk=zeros(N2-1,1);
  
  % Calculate expected periodogram for the dense th vector for all wavenumbers
  mnd=1;
  for ind=1:N2
      switch dth
          case 1
              th1=[thall{dth}(ind) th(2) th(3)];
          case 2
              th1=[th(1) thall{dth}(ind) th(3)];
          case 3
              th1=[th(1) th(2) thall{dth}(ind)];
      end
      Sks(ind,:)=maternos(k,th1,[],dim,1);
      if any(ind==thind)
          % if the current Matern parameter is one that we will plot profiles for,
          % take the opportunity to calculate the partial derivative and store the
          % periodogram separately
          dSk1s(mnd,:)=maternos(k,th1,dth,dim,1);
          Sk1s(mnd,:)=Sks(ind,:);
          mnd=mnd+1;
      else
      %Sks(ind-mnd+1,:)=maternos(k,th1,[],dim,1);
          dSk2s(ind-mnd+1,:)=maternos(ksub,th1,dth,dim,1);
          if ind>1
              r=Sks(ind,:)./Sks(ind-1,:)-1;
              %r=r(~isnan(r));r=r(~isinf(r));
              rSk(ind-1)=norm(r,2);
              rd=dSk2s(ind-mnd+1,:)./dSk2s(ind-mnd,:)-1;
              %rd=rd(~isnan(rd));rd=rd(~isinf(rd));
              rdSk(ind-1)=norm(rd,2);
          end
      end
  end
  
  % Let's normalize here
  Sks=(Sks'./S0)';
  switch dth
      case 1
          Sk1s=(Sk1s'./S0(thind))';
          dSk1s=(dSk1s'./S0(thind))';
      case 2
          Sk1s=(Sk1s'/S0)';
          dSk1s=(dSk1s'/S0)';
      case 3
          Sk1s=(Sk1s'./S0(thind))';
          dSk1s=(dSk1s'./S0(thind))';
  end
  Sk2s=(Sk2s'./S0d)';
  dSk2s=(dSk2s'./S0d)';
  
  % Parse dense Sks for wavenumbers of interest for the parameter dependent plots
  % and calculate the partial derivative
  for jnd=1:numcs
      temp=Sks(:,knd(jnd));
      Sk2s(:,jnd)=[temp(1:thind(1)-1); temp(thind(1)+1:thind(2)-1);...
                  temp(thind(2)+1:thind(3)-1); temp(thind(3)+1:end)];
                  %Sks(:,knd(jnd)); %maternos(ksub,th2,[],dim,1);
  end
  % Since we skip the thind indices when calculating the periodogram norm ratios,
  % lets set any zeros in the ratios to NaNs for plotting
  rSk(rSk==0)=NaN;
  rSk(find(isnan(rSk)==1)+1)=NaN;
  rdSk(rdSk==0)=NaN;
  
  % Calculate the central difference approximation for wavenumber dependence
  for ind=1:length(k)
      for jnd=2:N2-1
          % across columns, jnd: th(dth); down rows, ind: k
          % want central difference for changes in th(dth)
          cdiffSks(ind,jnd)=(Sks(jnd+1,ind)-Sks(jnd-1,ind))/(2*(thall{dth}(jnd+1)-thall{dth}(jnd-1)));
      end
  % Forward difference for first obs
      cdiffSks(ind,1)=(Sks(2,ind)-Sks(1,ind))/(thall{dth}(2)-thall{dth}(1));
  % Backward difference for last obs
      cdiffSks(ind,N2)=(Sks(N2,ind)-Sks(N2-1,ind))/(thall{dth}(N2)-thall{dth}(N2-1));
  end
  % We only need the central difference approximation for columns indexed
  % by thind and rows indexed by knd
  cdiffSk1s=cdiffSks(:,thind);
  cdiffSk2s=cdiffSks(knd,:);
  
  % The colors, marker sizes
  clrs=colormap(parula(4));
  ms=4;
  
  % Create the figure with three column plots
  figure(1);clf;
  [ah,ha,H]=krijetem(subnum(1,3));
  
  % Create the first suplot
  axes(ha(1))
  hold on
  for ind=1:numcs
      plot([ksub(ind) ksub(ind)],[-0.5 1.5],...
           ':', 'Color',[0.65 0.65 0.65],'LineWidth',0.1);
  end
  
  if ll
      pl1=semilogx(k,Sk1s,'-',...
                   k,dSk1s,'-.',...
                   k,2*cdiffSk1s,'o');
      pl1(end-2).MarkerSize=ms;
      pl1(end-1).MarkerSize=ms;
      pl1(end  ).MarkerSize=ms;
  else
      plot(k,Sk1s,'-',k,dSk1s,':',k,cdiffSk1s,'o');
  end
  xl(1)=xlabel('k, wavenumber (rad/km)');
  yl(1)=ylabel('normalized by S0');
  
  % Let's create a concise legend where the colors (values of parameter) and
  % marker style (function) only appear once each
  for ind=1:numcs
     % Specify the numcs different values by color
     LC1(ind)=plot(nan,nan,...
              'LineWidth',2,'DisplayName',...
              sprintf('$%s_%i$=%0.2f',labs{dth},ind,thsub{dth}(ind)));
  end
  % Specify the functions by line or marker style
  LC1(numcs+1)=plot(nan,nan,...
              'LineWidth',2,'Color','k','DisplayName',...
              sprintf('$S(k; %s_i)$',labs{dth}));
  LC1(numcs+2)=plot(nan,nan,...
              'LineStyle','-.','LineWidth',2,'Color','k','DisplayName',...
              sprintf('$\\partial S (k;%s_i) / \\partial %s$',...
                      labs{dth},labs{dth}));
  LC1(numcs+3)=plot(nan,nan,...
              'ko','DisplayName',...
              sprintf('$\\tilde{\\partial S}(k;%s_i) / \\partial %s$',...
                      labs{dth},labs{dth}));
  LCgnd1=legend(LC1,'Location','best','NumColumns',1,'Interpreter','latex');
  LCgnd1.ItemTokenSize=[12 12];
  
  % Clean up the axes
  yscl=0.025;
  xlim([k(1) max(k)]);
  ylim([min(min(dSk1s))-yscl max(max(Sk1s))+yscl])
  if ll
      xscale('log');
  end
  box on;
  set(gca,"ColorOrder",clrs(1:3,:))
  ytck=get(gca,'YTick');set(gca,'YTick',ytck(1:2:end));
  longticks;
  hold off
  
  % Create the middle subplot                                                     
  axes(ha(2))                                                                     
  hold on                                                                         
  for ind=1:numcs                                                                  
      vl2(ind)=plot([thsub{dth}(ind) thsub{dth}(ind)],[1e-10 1e5],...             
               ':','Color',[0.65 0.65 0.65],'LineWidth',0.5);                     
  end                                                                             
  if ll                                                                           
      pl2=loglog(thalld{dth},Sk2s,'-',...                                         
                 thalld{dth},dSk2s,'-.',...                                       
                 thall{dth}(2:end-1),2*cdiffSk2s(:,2:end-1),'o');                 
      pl2(end-2).MarkerSize=ms;                                                   
      pl2(end-1).MarkerSize=ms;                                                   
      pl2(end  ).MarkerSize=ms;                                                   
  else                                                                            
      plot(thall{dth},Sk2s,'-',...                                                
           thall{dth},dSk2s,'-.',...                                              
           thall{dth}(2:end-1),2*cdiffSk2s(2:end-1,:),'o');                       
  end                                                                             
                                                                                  
  % Let's create a concise legend where the colors (values of parameter) and      
  % marker style (function) only appear once each                                 
  for ind=1:numcs                                                                  
     % Specify the numcs different values by color                                 
     LC2(ind)=plot(nan,nan,...                                                    
              'LineWidth',2,'DisplayName',...                                     
              sprintf('$k_%i$=%0.1e',ind,ksub(ind)));                             
  end                                                                             
  % Specify the functions by marker/line style                                    
  LC2(numcs+1)=plot(nan,nan,...                                                    
              'LineWidth',2,'Color','k','DisplayName',...                         
              sprintf('$S(%s;k_i)$',labs{dth}));                                  
  LC2(numcs+2)=plot(nan,nan,...                                                    
              'LineWidth',2,'LineStyle','-.','Color','k','DisplayName',...        
              sprintf('$\\partial S (%s;k_i) / \\partial %s$',...                 
              labs{dth},labs{dth}));                                              
  LC2(numcs+3)=plot(nan,nan,...                                                    
              'ko','DisplayName',...                                              
              sprintf('$\\tilde{\\partial S}(%s;k_i) / \\partial %s$',...         
              labs{dth},labs{dth}));                                              
  LCgnd2=legend(LC2,'Location','best','NumColumns',1,'Interpreter','latex');      
  LCgnd2.ItemTokenSize=[12 12];                                                   
                                                                                  
  % Clean up the axes  
  xlabel(axlabs{dth})
  yscl=[2 1.1];
  ylim([min(min([dSk2s Sk2s])) max(max([dSk2s Sk2s]))].*yscl)
  if ll
      xscale('log');
      yscale('log');
  end
  box on;
  set(gca,"ColorOrder",clrs(1:3,:))
  ytck=get(gca,'YTick');set(gca,'YTick',ytck(1:2:end));
  longticks;
  hold off
  
  % Make a title with information relevant to creating the periodograms
  switch dth
      case 1
          tistr=sprintf('$\\nu=%0.1g,\\rho=%0.1g$',th(2),th(3));
      case 2
          tistr=sprintf('$\\sigma^2 = %0.1g, \\rho = %0.1g$',th(1),th(3));
      case 3
          tistr=sprintf('$\\sigma^2 = %0.1g, \\nu = %0.1g$',th(1),th(2));
  end
  ti=title(sprintf('%s;  %ix%i grid;  %ix%i km',...
           tistr,p.NyNx,p.NyNx.*p.dydx/1e3),'Interpreter','latex');
  movev(ti,0.01);
  
  % Make the final subplot
  axes(ha(3))
  hold on
  for ind=1:numcs
      vl3(ind)=plot([thsub{dth}(ind) thsub{dth}(ind)],[1e-20 max(max([rSk; rdSk]))+1],...
               ':','Color',[0.65 0.65 0.65],'LineWidth',0.1);
  end
  set(gca,"ColorOrder",...
      [0 0.4470 0.7410;0.9290 0.6940 0.1250;0.6350 0.0780 0.180])
  if ll
      pl3=loglog(thall{dth}(2:end),rSk,'-',...
                   thall{dth}(2:end),rdSk,'-.');
  else
      plot(thall{dth}(2:end),rSk,'-',...
           thall{dth}(2:end),rdSk,'-.');
  end
  
  % Let's create a nice legend
  LC3(1)=plot(nan,nan,...
              'LineWidth',2,'DisplayName',...
              sprintf('$||\\frac{S(k; %s_{i+1})}{S(k; %s_{i})} - 1||_2$',...
                      labs{dth},labs{dth}));
  LC3(2)=plot(nan,nan,...
              '-.','LineWidth',2,'DisplayName',...
     sprintf('$||\\frac{ \\partial S(k; %s_{i+1})/ \\partial %s}{ \\partial S(k;%s{i} / \\partial %s } - 1 ||_2$',...
               labs{dth},labs{dth},labs{dth},labs{dth}));
  LCgnd3=legend(LC3,'NumColumns',1,'Interpreter','latex');
  LCgnd3.ItemTokenSize=[12 12];
  
  % Clean up the axes
  xlabel(axlabs{dth})
  yscl=[0.9 1.1];
  ylim([min(min([rSk rdSk])) max(max([rSk rdSk]))].*yscl)
  if ll
      xscale('log');
      yscale('log');
  end
  box on;
  set(gca,"ColorOrder",clrs(1,:))
  longticks;
  hold off
  
  % Save the figure
  switch dth
      case 1
          tith=[th(2) th(3)];
      case 2
          tith=[th(1) th(3)];
      case 3
          tith=[th(1) th(2)];
  end
  figdisp([],sprintf('demo1_%3.3i_%3.3i_%idth_%i_%i',p.NyNx,dth,tith),[],1);
end
