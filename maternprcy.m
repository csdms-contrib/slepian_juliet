function varargout=maternprcy(theta,prc,meth,xver)
% [Yprc,yrange,Cy,fzini,aneval0toY,trapeval0toY,yfun]=MATERNPRCY(theta,prc,meth,xver)
%
% Given a Matern parameter vector theta [s2 nu rho], returns the
% lag distance Y for which $int_\0^(Y) y C(y;\theta) dy$
% equals a target percentage of $int_\0^\inf C(y;\theta) dy$. 
% We can additionally calculate these integrals for C(y;\theta) alone.
%
% INPUT
% 
% theta   vector of Matern parameters [s2 nu rho]
% prc     target percentile of overall area under Cy [default: 50]
% meth    method for selecting initial guess for root-finding, 1 is quick, 2 is
%         better
% xver    0 no extra verification [default]
%         1 extra verification via evaluation of MATERNOSY with the trapezoidal
%           rule over a discretized range of lag-distances
% 
% OUTPUT
%
% Yprc         lag-distance corresponding to the target percentage, PRC
% yrange       some trial lag distance range, used to create fzini
% Cy           the actual Matern function
% fzini        the initial guess used for root-finding, selected based on meth
% aneval0toY   analytic function evaluation for the area under the curve
% trapeval0toY numerical integration over yrange up to Yprc, only if xver==1
%
% EXAMPLE
%
% [Yprc,yrange,Cy,fzini]=maternprcy([20  2.5 35000],50,2,1); 
% 
% maternprcy('demo1') % Evaluation, with and without verification
% maternprcy('demo2') % A plot of the evaluation
% maternprcy('demo3') % Calculations over a sequence of percentages 
% maternprcy('demo4') % Calculations over a sequence of smoothness parameters 
%
% SEE ALSO
%
% MATERNOSY, MATERNPRC, SYM, SYMS, FZERO, TRAPEZE
%
% TESTED R2017a (9.2.0.538062) 64-bit (glnxa64), 
%        R2024b Update 3 (24.2.0.2806996) 64-bit (glnxa64)
%
% Last modified by olwalbert-at-princeton.edu, 01/16/2026

% Default values
defval('prc',50)
defval('meth',2)
defval('xver',0)
defval('theta',[10 0.75 250])
% Still need one for the recursive part where xver is set to 0
defval('yrange',NaN);
defval('Cy',NaN);
defval('fzini',[NaN NaN]);

% Not a demo!
if ~isstr(theta)
    % shunt the 0 and 100 percentiles but note
    % that in those cases you only ever get one output
    if prc==0
      Yprc=0; 
      varns={Yprc,yrange,Cy,fzini};
      varargout=varns(1:nargout);
      return
    elseif prc==100
      % For Cy,
      % Yprc=theta(1)*pi*theta(3)/2/sqrt(theta(2))*gamma(1/2+theta(2))/gamma(theta(2));
      % For yCy,
      Yprc=s2/2*(pee*rho)^2;
      varns={Yprc,yrange,Cy,fzini};
      varargout=varns(1:nargout);
      return
    end

    % If you'll be needing a yrange output or input, make it here
    if nargout>1 || xver==1
      % Start with a trial lag-distance range 
      if theta(3)>100
        yrange=linspace(0,max([10^(ceil(log10(theta(3)))+1) 5e3]),1e4);
      elseif theta(3)>10
        yrange=linspace(0,max([10^(ceil(log10(theta(3)))+1) 5e2]),1e4);
      else
        yrange=linspace(0,max([10^(ceil(log10(theta(3)))+1) 5e1]),1e4);
      end
    end

    % And this is the theoretical y-axis intercept
    C0=theta(1);

    % The symbolic expressions
    % MATLAB pre-2018a is needed for the following syntax to run
    % Cys=sym(...
    %     's2*2^(1-nu)/gamma(nu)*(2*nu^(1/2)*y/(pee*rho))^nu*besselk(nu,2*nu^(1/2)*y/(pee*rho))');
    syms s2 nu rho positive
    syms y Y pee prcnt positive
    Cys=s2*2^(1-nu)/gamma(nu)*(2*nu^(1/2)*y/(pee*rho))^nu*besselk(nu,2*nu^(1/2)*y/(pee*rho));

    % MATLAB cannot evaluate integrals for the bessel function. If it could,
    % we would have done verifications involving the following
   	% CintY=simplify(int(Cys,'y',0,Y));
    % Cintinf=simplify(int(Cys,'y',0,Inf));
    % This will fail, too
    % Cym=@(yd) double(subs(Cys,[y,s2,nu,rho,pee],[yd,theta(1),theta(2),theta(3),pi]));
    % Cintinf=integral(Cym,0,Inf)
   	% and then explicitly check by subbing things like
    % % difer(subs(Cintinf-auc0toInf),[],[],NaN)
    % % difer(subs(Cintinf-auc0toInf,[Y s2 nu rho pee],[rand theta pi]),[],[],NaN)
    % % difer(subs(simplify(Cintinf-auc0toInf),[Y s2 nu rho pee],[rand theta pi]),[],[],NaN)
        
    % The area under the curve is the ONE-dimensional integral in the
    % isotropic case
    % A = \int_0^Y y C(y) dy
    % Keep the variable notation for symmetry with MATERNPRC, but do notice 
    % that we do not have to worry about dimension for the isotropic spatial 
    % covariance
        
    % The total area under the curve from 0 to inf for Cy alone is found by 
    % applying G&R 6.561.16
    % auc0toInf=s2*pee^(3/2)*rho/(2*nu^(1/2))*gamma(1/2+nu)/gamma(nu);
    % For yCy, this is simply
    auc0toInf=s2/2*(pee*rho)^2;

    % The area under the curve from 0 to Y will be between 0 and auc0toInf
    arguh=nu^(1/2)*Y/(pee*rho);
    % For Cy alone, this was complicated and is documented on Overleaf
    % auc0toY=s2*(nu^2+nu)*csc(pee*nu)/2^(nu)*(2*arguh)^(2*nu+2)*...
    %         (gamma(nu+1)/gamma(-nu+1)*hypergeom([nu+1],[-nu+1,nu+2],arguh^2)+...
    %         hypergeom([],[nu+2],arguh^2));
    % For yCy, applying IBP or eq. (7.14.1.3) of Erdelyi 1953 gives
    auc0toY=1/2*pee*rho*s2*(pee*rho-2*(nu^(1/2*(-1+nu)))*pee^(-nu)*Y*(Y/rho)^nu*...
            besselk(1+nu,2*arguh)/gamma(nu));

    % Recall that to substitute in values for the symbolic parameters and 
    % evaluate the expression, we use the following
    % double(subs(auc0toY,[Y s2 nu rho pee],[1e8 theta pi]));

    % In contrast to MATERNPRC, we cannot solve for Y, so we have the option of
    % using numerical bisection of the analytical system
    % prcnt*auc0toInf-auc0toY
    yfun=prcnt/100*auc0toInf-auc0toY;
    % that we will verify by summing with TRAPEZE

    % Find the root(s) of yfun using a bisection algorithm. First, set up the 
    % symbolic function to evaluate and solve
    % Note that auc0toY is not defined for integer values of nu, so we will have
    % to add a small perturbation
    if mod(theta(2),1)~=0
      % nu is not an integer
      Y2solvemf=@(Yd,s2d,nud,rhod,peed,prcntd) ...
         real(double(subs(yfun,[Y,s2,nu,rho,pee,prcnt],[Yd,s2d,nud,rhod,peed,prcntd])));
    else
      % nu is an integer; sometimes we end up with complex values
      Y2solvemf=@(Yd,s2d,nud,rhod,peed,prcntd) ...
        real(double(subs(yfun,[Y,s2,nu,rho,pee,prcnt],[Yd,s2d,nud+1e-3,rhod,peed,prcntd]))+...
        double(subs(yfun,[Y,s2,nu,rho,pee,prcnt],[Yd,s2d,nud-1e-3,rhod,peed,prcntd]))/2);
    end

    % Then assign a proper guess relative to the actual spatial covariance
    Cy=maternosy(yrange,theta);
    % Check that the intercept is correct to 1e-7 or thereabouts
    difer(maternosy(0,theta)-C0,7,[],NaN);
    % Find the halfway point of the curve to use as a starting guess for Y
    Cygues=yrange(find(diff(sign(Cy/C0-0.5))));  
    % Or, find the largest lag-distance where the normalized Cy curve
    % reaches a low value (to being the interesting nook) 0.01
    try
      Cygues=yguess(yrange,Cy,prc,meth,yfun,theta);
    catch
      % If this failed, try to automate the troubleshooting
      try
        if meth==1
          meth=2;
          Cygues=yguess(yrange,Cy,prc,meth,yfun,theta);
        else
          yrange2=linspace(1,max([10^(ceil(log10(theta(3)))+1) 5e3]),yrange(end)*100);
          Cygues=yguess(yrange2,Cy,prc,meth,yfun,theta);
        end
      catch
        % and if that failed, intervene
        keyboard
      end
    end
    switch meth
      case 1
        fzini=Cygues(1);
      case 2
        if size(Cygues,1)>1
          % Make sure that the bracket is set up correctly when we encounter
          % multiple roots
          keyboard
        end
        fzini=Cygues;
    end
    % Options for exploring the root-finder
    options = optimset('PlotFcns','optimplotfval','Display','iter');
    options = [];
    try
      [Yprc,Yval,fzexit,fzout]=fzero(@(Y)...
               Y2solvemf(Y,theta(1),theta(2),theta(3),pi,prc),fzini,options);
    catch
      try
        % a common reason to end up here is that the fzini bracket needs to be
        % larger, so give it some room to breathe
        fzini=[fzini(1)*0.8 fzini(2)*1.2];
        [Yprc,Yval,fzexit,fzout]=fzero(@(Y)...
                 Y2solvemf(Y,theta(1),theta(2),theta(3),pi,prc),fzini,options);
      catch
        % intervene 
        keyboard
      end
    end
    Yprc=abs(Yprc);
    
    % More output requested, but you hadn't needed them: make them 
    if nargout>2 && [exist('Cy','var')==0 || strcmp(class(Cy),'sym') || any(isnan(Cy))]
      Cy=maternosy(yrange,theta);
    end
    if nargout>3 && [exist('Cygues','var')==0 || any(isnan(Cygues))]
      Cygues=yguess(yrange,Cy,prc,1,yfun,theta);
    end
        
    % If extra verification is requested, compare the calculated results
    if xver == 1
        % Check via explicit numerical integration on the yrange
        trapeval0toInf=trapeze(yrange(:),yrange(:).*maternosy(yrange(:),theta));
        % equivalent built-in function
        % Q=trapz(yrange(:),maternosy(yrange(:),theta));
        aneval0toInf=double(subs(auc0toInf,[s2 nu rho pee],[theta pi]));
        disp(sprintf('\n0 to "Inf" - %g, TRAPEZE %g Relative Difference %g',...
		     aneval0toInf,trapeval0toInf,...
             (aneval0toInf-trapeval0toInf)/trapeval0toInf))

        % Check "Yprc": Compare the calculated area under the curve from 0 to Y
        % to the area that Y corresponds to as a percentage of the total area
        ytoY=yrange(yrange(:)<=Yprc); ytoY=ytoY(:);
        trapeval0toY=trapeze(ytoY,ytoY(:).*maternosy(ytoY,theta));
        if mod(theta(2),1)~=0
          % nu is not an integer
          aneval0toY=double(subs(auc0toY,[Y s2 nu rho pee],[Yprc theta pi]));
        else
          % nu is an integer and we have to average over perturbed values
          aneval0toY=(double(subs(auc0toY,...
                       [Y s2 nu rho pee],[Yprc theta(1) theta(2)-1e-3 theta(3) pi]))+...
                       double(subs(auc0toY,...
                       [Y s2 nu rho pee],[Yprc theta(1) theta(2)+1e-3 theta(3) pi])))/2;
        end
        disp(sprintf('\n0 to Y ( %g ) for %g%s - %g, TRAPEZE %g Relative Difference %g\n',...
		     Yprc,prc,'%',aneval0toY,trapeval0toY,...
             (aneval0toY-trapeval0toY)/trapeval0toY))
        disp(sprintf('\n analytical %g%%, TRAPEZE %g%%, goal %g%%\n',...
             aneval0toY/aneval0toInf*100,trapeval0toY/trapeval0toInf*100,prc))
    else
      aneval0toY=[];
      trapeval0toY=[];
      aneval0toInf=[];
      trapeval0toInf=[];
    end
    % Provide as much output as requested
    varns={Yprc,yrange,Cy,fzini,aneval0toY,trapeval0toY,aneval0toInf,trapeval0toInf,yfun};
    varargout=varns(1:nargout);
  elseif strcmp(theta,'demo1')
    % with and without verification; not so interesting
    prc=rand*100;
    [Yprc,yrange,Cy,fzini]=maternprcy([],prc,1,0); disp(sprintf('%g',Yprc))
    [Yprc,yrange,Cy,fzini]=maternprcy([],prc,1,1); disp(sprintf('%g',Yprc))
    [Yprc,yrange,Cy,fzini]=maternprcy([],prc,2,0); disp(sprintf('%g',Yprc))
    [Yprc,yrange,Cy,fzini]=maternprcy([],prc,2,1); disp(sprintf('%g',Yprc))
  elseif strcmp(theta,'demo2')
    % make a quick plot of the normalized spatial covariance with its
    % PRC percentile marked
    prc=rand*100; xver=1;%round(rand);
    meth=2;
    thet=[1 max(round(rand*10,1),0.5) max(round(rand*10000),0.00001)];
    disp(sprintf('s2 %g nu %g rho %g',thet))
    disp(sprintf('Trying the %gth percentile',prc))
    [a,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver);
    plotcy(a,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun)
  elseif strcmp(theta,'demo3')
    thet=[1 rand(1)*10 max(round(rand*10000),0.00001)];
    thet=[1 0.5 250];
    disp(sprintf('s2 %g nu %g rho %g',thet))
    xver=1;
    meth=2;
    %for prc=1:3:99
    for prc=5:10:100
      [a,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver); 
      plotcy(a,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun)
      pause(0.5)
    end
  elseif strcmp(theta,'demo4')
    prc=50;
    xver=1;
    meth=2;
    ps=0.5;
    disp(sprintf('Isotropic Matern %ith percentile xver %i',...
  	       round(prc),xver))
    thet=[1 0.33 250];
    [Yprc,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver); 
    disp(sprintf('%g',Yprc));
    plotcy(Yprc,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun); 
    pause(ps)
  
    thet=[1 0.5 250];
    [Yprc,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver); 
    disp(sprintf('%g',Yprc));
    plotcy(Yprc,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun);
    pause(ps)

    thet=[1 1 250];
    [Yprc,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver); 
    disp(sprintf('%g',Yprc));
    plotcy(Yprc,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun);
    pause(ps)
  
    thet=[1 1.5 250];
    [Yprc,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver); 
    disp(sprintf('%g',Yprc));
    plotcy(Yprc,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun);
    pause(ps)
    
    thet=[1 2.5 250];
    [Yprc,b,c,d,e,f,g,h,yfun]=maternprcy(thet,prc,meth,xver); 
    disp(sprintf('%g',Yprc));
    plotcy(Yprc,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun);
    pause(ps)
  end
end

function plotcy(a,b,c,d,e,f,g,h,thet,prc,meth,xver,yfun)
  % a: Yprc, b: yrange, c: Cy, d: cygues, e: aneval0toY, f: trapeval0toY,
  % g:aneval0toInf, h: trapeval0toInf
  sv=1;
  figure(2);clf
  ah(1)=subplot(121);
  % plot the spatial covariance
  Cyp(1)=plot(b,c,'LineW',2,'DisplayName','C(y)');
  hold on
  Cyp(2)=plot(b(b<=a),c(b<=a),'r','LineW',2,'HandleVisibility','off');
  % plot the lag-distance corresponding to prc
  plot([a a],ylim,'r-','LineW',1,'DisplayName',sprintf('%5.2f %%',prc))
  offs=0.05;
  [~,lgc]=min(abs(b-a));
  t(1)=text(a+offs*10^(log10(a)+1),0.5*c(lgc)+offs,...
       sprintf('auc an.: %5.2f\nauc num.: %5.2f\n%s=%5.2f',e,f,'y_\alpha',a),...
       'Interpreter','tex');
  grid on
  if xver~=1
    plot(d(1),d(2),'o','MarkerF','r','MarkerE','r','DisplayName','guess')
    hold off
  end
  legend
  ylim([-0.01 1.01])
  xlabel('lag-distance y')
  ylabel('C(y)')
  % now evaluate the symbolic function we solved for along the decimated yrange
  syms s2 nu rho y Y pee prcnt positive
  decm=1e2;
  for ind=1:floor(numel(b)/decm);
    yfunev(ind)=double(subs(yfun,[Y s2 nu rho pee prcnt],[b(decm*(ind-1)+1) thet pi prc]));
  end
  ah(2)=subplot(122);
  % Plot the evaluated function to see its intercept 0 at the solved
  % lag-distance
  plot(b(1:decm:end),yfunev);ylim([-max(yfunev) max(yfunev)])
  syfunev=sign(yfunev);
  schk=syfunev(1:end-1).*syfunev(2:end);
  snd=find(schk<0);
  disp(sprintf('There are %i roots near %5.2f %5.2f',sum(schk<0),b(decm*(snd-1)+1)))
  hold on
  plot([a a],ylim,'r-','LineW',1,'DisplayName',sprintf('%5.2f %%',prc))
  t(2)=text(a+offs*10^(log10(a)),offs*10^floor(log10(max(yfunev))+1),...
       sprintf('%s=%5.2f','y_\alpha',a),'Interpreter','tex');
  yline(0)
  grid on
  xlabel('lag-distance y')
  ylabel('\alpha \int_0^\infty C(r) dr - \int_0^Y C(r) dr')
  sgtitle(sprintf('%5.2fth percentile for Matern [%5.2f %5.2f %5.2f], meth %i, xver %i',...
        prc,thet,meth,xver))
  if sv
    saveas(gcf,...                                                               
           sprintf('maternprcy_plotcy_s2%i_10nu%i_rh%i_%i.eps',...           
           round(thet.*[1 10 1]),round(prc)),'epsc') 
  end
end

function yg=yguess(yr,Cy,prc,meth,yfun,thet)
  defval('meth',1)
  switch meth
    case 1
      % This is consistent with KGUESS from MATERNPRC which finds a point near
      % the inflection of the spatial covariance as the initialization. It is 
      % faster but has worse performance than case 2.
      % There is art in picking this level - if FZERO complains, change this.
      thepoint=0.01;
      % Here's a trick 
      cands=yr(find(diff(...
          sign(Cy/max(Cy)-thepoint))));
      if prc<50
        yg(1)=min(cands);
      else
        yg(1)=max(cands);
      end
      yg(2)=thepoint;
    case 2
      % Evaluate the symbolic function for a decimated range of lag-distances to
      % approximate the root and provide an interval for the initizialization of
      % FZERO.
      syms s2 nu rho y Y pee prcnt positive
      decm=1e2;
      for ind=1:floor(numel(yr)/decm);
        yfunev(ind)=double(subs(yfun,[Y s2 nu rho pee prcnt],[yr(decm*(ind-1)+1) thet pi prc]));
      end
      syfunev=sign(yfunev);                                                         
      schk=syfunev(1:end-1).*syfunev(2:end);                                        
      snd=find(schk<0);   
      if sum(schk<0)<1
        % No root was found. This is likely due to a too short or too coarse yr.
        keyboard
      end
      % To bracket the interval(s) that FZERO will search, we must set the
      % bounds so that they evaluate to opposite signs. We will provide values 
      % that are adequately spaced around the approximate root.
      for jnd=1:numel(snd)
        yg(jnd,:)=[yr(decm*(snd(jnd)-1)) yr(decm*(snd(jnd)))];
      end
  end
end
