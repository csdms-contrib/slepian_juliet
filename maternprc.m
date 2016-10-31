function varargout=maternprc(theta,prc,meth,xver)
% [Kprc,krange,Sk,kSkgues]=MATERNPRC(theta,prc,meth,xver)
%
% Given a two-dimensional Matern parameter vector theta [s2 nu rho] 
% returns the wavenumber K for which $int_\0^K S(k;\theta) k dk$
% equals a target percentage of $int_\0^\inf S(k;\theta) k dk$
%
% INPUT:
% 
% theta   vector of Matern parameters [s2 nu rho]
% prc     target percentile of overall area under Sk [default: 50]
% meth    0 symbolically build up to numerical/analytic methods (slow!) 
%         1 numerical (using an initial guess)
%         2 analytical [default, fast, and preferred]
% xver    0 no extra verification [default]
%         1 extra verification via the symbolic method, and after the
%           fact for the other methods
% 
% OUTPUT:
%
% Kprc       wavenumber corresponding to target percentage
% krange     some trial wavenumber range, used to create kSkgues
% Sk         the actual Matern function
% kSkgues    the half-Sk wavenumber, used as a guess for meth==1
%
% NOTE: The FZERO function used requires an initial guess, picked from a
%       logarithmically spaced set of wavenumbers that should be good for
%       most Matern Sk curves, but this is not guaranteed.
%
% EXAMPLE:
%
% maternprc('demo1') % Just a bunch of default evaluations
% maternprc('demo2') % A bunch of random evaluations with some plots
% maternprc('demo3') % A bunch of sequential evaluations
% maternprc('demo4') % Gabe no-like values
%
% TO DO: shunt the 0th and 100th 
%
% SEE ALSO:  MATERNOS, MATERNOSP, MATERNOSY, FZERO, SYM, SYMS, TRAPEZE
%
% Last modified by gleggers-at-princeton.edu, 03/19/2014
% Last modified by fjsimons-at-princeton.edu, 10/31/2016

% Default values
defval('prc',50)
defval('meth',2)
defval('xver',0)
defval('theta',[0.022  2.5 35000])
% Still need one for the recursive part where xver is set to 0
defval('krange',NaN);
defval('Sk',NaN);
defval('kSkgues',[NaN NaN]);

% Not a demo!
if ~isstr(theta)
    % Regardless of the method, shunt the 0 and 100 percentiles but note
    % that in those cases you only ever get one output
    if prc==0
      Kprc=0; 
      varns={Kprc,krange,Sk,kSkgues};
      varargout=varns(1:nargout);
      return
    elseif prc==100
      Kprc=theta(1);
      varns={Kprc,krange,Sk,kSkgues};
      varargout=varns(1:nargout);
      return
    end

    % If you'll be needing a krange output or input, make it here
    if nargout>1 || meth==1 || xver==1
      % Start with a trial wavenumber range 
      krange=logspace(-9,-1,1e4);
    end

    % And this is the theoretical y-axis intercept
    S0=[theta(1)*theta(3)^2/4*pi];

    % Symbolic math method, and extra verification
    if meth==0 || xver==1
        % Matern isotropic 2D spectral density, Simons & Olhede (2013), eq. (72)
        syms s2 nu rho positive
	syms k K pee prcnt positive
        Sk=sym(...
	    's2*nu^(nu+1)*4^nu/pee/(pee*rho)^(2*nu)*(4*nu/(pee*rho)^2+k^2)^(-nu-1)');
        
        % The area under the curve is the TWO-dimensional integral that
        % in polar coordinates is a ONE-dimensional integral 
        % A = \int_0^K \int_0^K S(k) d\bk
        %   = 2\pi \int_0^K S(k) k dk
        % The integrand of the equivalent ONE-dimensional integral
        Skk=2*pee*Sk*k;
        
        % The total area under the curve from 0 to inf simply IS s2...
        Sintinf=simplify(int(Skk,'k',0,Inf));
	
	% ... which we should explicitly check by subbing anything
	difer(subs(Sintinf-s2),[],[],NaN)
        
        % The area under the curve from 0 to K ranges between 0 and s2
	SintK=simplify(int(Skk,'k',0,K));

	% This is the evaluated Sintk spelled out...
	spl1=2*4^nu*nu^(nu + 1)*s2/(pee*rho)^(2*nu)*...
	     [((pee^2*rho^2)/(4*nu))^nu/(2*nu)-1/(2*nu*(K^2+(4*nu)/(pee^2*rho^2))^nu)];
	% ... which we should explicitly check by subbing anything (even
        % nothing!) into what should be exactly, symbolically zero
	difer(subs(SintK-spl1,symvar(spl1)),[],[],NaN)
		
	% And this is how we simplify it even further, by hand
	spl2=s2*[1-2^(2*nu)*nu^nu/(K^2*pee^2*rho^2 + 4*nu)^nu];
	% Divide by s2 since we'll have a right hand side containing s2
	spl3=   [1-2^(2*nu)*nu^nu/(K^2*pee^2*rho^2 + 4*nu)^nu];
	
	% you can verify that these two representations are identical by
        % substitution, e.g.  
	% ... this should be exactly, symbolically zero
	difer(subs(symvar(spl2)-symvar(spl1)),[],[],NaN)
	% ... this will be exactly, numerically zero (I need to evaluate)
        % since Matlab isn't smart enough to realize they are actually,
        % symbolically zero, except perhaps after simplification
	difer(subs(spl1-spl2,[K s2 nu rho pee],[rand theta pi]),[],[],NaN)
	% ... turns out that simplify brings out the symbolic zero
	difer(subs(simplify(spl1-spl2)),[],[],NaN)

	% or you can turn them into symbolic functions, e.g.
	spl1f=symfun(spl1,[K s2 nu rho pee]);
	spl2f=symfun(spl2,[K s2 nu rho pee]);
	% Watch the order of the input variables - symvar lists
        % alphabetically
	rando=rand;
	difer(eval(spl1f(rando,theta(1),theta(2),theta(3),pi))-...
	      eval(spl2f(rando,theta(1),theta(2),theta(3),pi)),...
	      [],[],NaN)

	% Solve for a certain percentile of s2; divide by s2 then solve this
	K2solv=symfun(spl3-prcnt,[K prcnt nu rho pee]);
	
        % We can write the solution out by hand - RB X page 25
        Ksolved=symfun(2/(pee*rho)*sqrt(nu*[(1-prcnt)^(-1/nu)-1]),[prcnt nu rho pee]);
    end

    % Now cycle through the methods
    switch meth
     case 0
      % Symbolic derivation of the numerical and the analytical method
      % If you made it all the way here, you can just evaluate the
      % analytical solution...
      Kprc=eval(Ksolved(prc/100,theta(2),theta(3),pi));
      % Verfiy that you have just found the solution by plugging in
      difer(eval(K2solv(Kprc,prc/100,theta(2),theta(3),pi)),[],[],NaN)
      % ... or you can also solve the symbolic expression numerically
      % after converting it to an anonymous function handle
      K2solvemf=matlabFunction(K2solv);
      % Then solve numerically using FZERO to be sure
      % in this case we simply take Kprc as the initial guess...
      Kprc2=fzero(@(K) K2solvemf(K,theta(2),pi,prc/100,theta(3)),Kprc);
      % Check that they're both consistent
      difer(Kprc-Kprc2,[],[],NaN)
     case 1
      % Direct route to the numerical method...
      K2solvemf=@(K,nu,pee,prcnt,rho) ...
		-prcnt-2.0.^(nu.*2.0).*nu.^nu.*(nu.*4.0+K.^2.*pee.^2.*rho.^2).^(-nu)+1.0;

      % ... but now need proper guess
      % Here is the actual Matern function 
      Sk=maternos(krange,theta);
      % Check that the intercept is correct to 1e-7 or thereabouts
      difer(maternos(0,theta)-S0,7,[],NaN);
      % Find the halfway point of the curve to use as a starting guess K
      % kSkgues=krange(find(diff(sign(Sk/S0-0.5))));  
      % Find the largest wavenumber where the normalized kSk curve
      % reaches a low value (to being the interesting nook) 0.01
      kSkgues=kguess(krange,Sk,prc);

      % What is the solution in the end? One could force it all to be
      % positive - in earlier experiments we found the right number, only
      % the wrong sign (we doing other things differently...)
      Kprc=abs(fzero(@(K) K2solvemf(K,theta(2),pi,prc/100,theta(3)),kSkgues(1)));
     case 2
      % Direct route to the analytical method
      Kprc=2/(pi*theta(3))*sqrt(theta(2)*[(1-prc/100)^(-1/theta(2))-1]);
    end
    
    % More output requested, but you hadn't needed them: make them 
    if nargout>2 && [exist('Sk','var')==0 || strcmp(class(Sk),'sym') || any(isnan(Sk))]
      Sk=maternos(krange,theta);
    end
    if nargout>3 && [exist('kSkgues','var')==0 || any(isnan(kSkgues))]
      % kSkgues=krange(find(diff(sign(Sk/S0-0.5))));
      kSkgues=kguess(krange,Sk,prc);
    end
        
    % If extra verification is requested, compare the calculated results
    if xver == 1
        % Check "Sintinf" via explicit numerical integration on the krange
        cSintinf=trapeze(krange(:),2*pi*maternos(krange(:),theta).*krange(:));
        disp(sprintf('\nS2 %g TRAPEZE %g Difference %g',...
		     theta(1),cSintinf,theta(1)-cSintinf))

        if meth~=0
	  % Check "Kprc"via the "other" method (of course meth==0 is OK)
	  Kprc2=maternprc(theta,prc,meth+(-1)^(meth+1),0);
	  disp(sprintf('METHOD %g OTHER METHOD  %g Difference %g',...
		       Kprc,Kprc2,Kprc-Kprc2))
	end
	  
        % Check "Kprc": Compare the calculated area under the curve from 0 to K
        % to the area that K corresponds to as a percentage of the total area
        ktoK=krange(krange(:)<=Kprc); ktoK=ktoK(:);
        cSintK=trapeze(ktoK,2*pi*maternos(ktoK,theta).*ktoK);
        disp(sprintf('%s0K %g TRAPEZE %g Difference %g\n',...
		     '%',prc/100*theta(1),cSintK,prc/100*theta(1)-cSintK))
    end

    % Provide as much output as requested
    varns={Kprc,krange,Sk,kSkgues};
    varargout=varns(1:nargout);
    
elseif strcmp(theta,'demo1')
  prc=rand*100;
  [Kprc,krange,Sk,kguess]=maternprc([],prc,0,1); disp(sprintf('%g',Kprc))
  [Kprc,krange,Sk,kguess]=maternprc([],prc,1,1); disp(sprintf('%g',Kprc))
  [Kprc,krange,Sk,kguess]=maternprc([],prc,2,1); disp(sprintf('%g',Kprc))

  [Kprc,krange,Sk,kguess]=maternprc([],prc,0,0); disp(sprintf('%g',Kprc))
  [Kprc,krange,Sk,kguess]=maternprc([],prc,1,0); disp(sprintf('%g',Kprc))
  [Kprc,krange,Sk,kguess]=maternprc([],prc,2,0); disp(sprintf('%g',Kprc))
  
elseif strcmp(theta,'demo2')
  prc=rand*100; meth=randi(3)-1; xver=round(rand);
  thet=[max(round(rand*1000)/1000,0.0001) max(round(rand*10),0.5) max(round(rand*10000),0.00001)];
  disp(sprintf('s2 %g nu %g rho %g',thet))
  disp(sprintf('Trying the %gth percentile using method %i',prc,meth))
  [a,b,c,d]=maternprc(thet,prc,meth,xver);
  plotsk(a,b,c,d,prc,meth,xver)
elseif strcmp(theta,'demo3')
  thet=[max(round(rand*1000)/1000,0.0001) max(round(rand*10),0.5) max(round(rand*10000),0.00001)];
  disp(sprintf('s2 %g nu %g rho %g',thet))
  meth=randi(3)-1; xver=round(rand);

  for prc=1:3:99
    [a,b,c,d]=maternprc(thet,prc,meth,xver); 
    plotsk(a,b,c,d,prc,meth,xver)
    pause
  end
elseif strcmp(theta,'demo4')
  prc=rand*100; meth=randi(3)-1; xver=round(rand);
  
  prc=50;
  
  disp(sprintf('2-D Matern %ith percentile method %i xver %i',...
	       round(prc),meth,xver))
  
  [Kprc,b,c,d]=maternprc([0.1636 0.0049 8.3618e+03],prc,meth,xver); 
  disp(sprintf('%g',Kprc)); plotsk(Kprc,b,c,d,prc,meth,xver); pause

  [Kprc,b,c,d]=maternprc([0.0160 0.1146 5.9485e+04],prc,meth,xver); 
  disp(sprintf('%g',Kprc)); plotsk(Kprc,b,c,d,prc,meth,xver); pause

  [Kprc,b,c,d]=maternprc([0.0046 0.2629 3.6451e+04],prc,meth,xver); 
  disp(sprintf('%g',Kprc)); plotsk(Kprc,b,c,d,prc,meth,xver); pause
  
  [Kprc,b,c,d]=maternprc([7.5641e-04 0.3199 1.5281e+05],prc,meth,xver); 
  disp(sprintf('%g',Kprc)); plotsk(Kprc,b,c,d,prc,meth,xver); pause
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotsk(a,b,c,d,prc,meth,xver)
clf
Skp=semilogx(b,c/max(c),'LineW',2); hold on
Skkp=semilogx(b,c(:).*b(:)/max(c(:).*b(:)),'LineW',2); hold on
set(Skkp,'Color','r')
plot([a a],ylim,'r-','LineW',1); grid on
if meth~=0 && [meth~=2 && xver~=1]
  plot(d(1),d(2),'o','MarkerF','r','MarkerE','r'); hold off
end
title(sprintf('2-D Matern %5.2fth percentile method %i xver %i',prc,meth,xver))
ylim([-0.1 1.1])
xlabel('wavenumber k')
ylabel('S(k)/S(0) and normalized S(k)k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kg=kguess(kr,Sk,prc)

% There is art in picking this level - if FZERO complains, change this
thepoint=0.01;

% Here's a trick 
cands=kr(find(diff(...
    sign(kr(:).*Sk/max(kr(:).*Sk)-thepoint))));
if prc<50
  kg(1)=min(cands);
else
  kg(1)=max(cands);
end

kg(2)=thepoint;
