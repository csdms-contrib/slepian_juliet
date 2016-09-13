function eggers1
% EGGERS1
%
% The influence of unknown COVARIANCE on the brute-force estimation of
% the VARIANCE of a realization of the isotropic Matern process.
% You must open up MATLABPOOL to do multiple simulations per processor!
%
% SEE ALSO:
%
% SIMULOSL('demo4')
%
% Last modified by fjsimons-at-alum.mit.edu 05/08/2016

% Metric conversion
mfromkm=1000;

% Matern parameters, in standard units
% Variance of the "topography" [m^2]
s2=[1^2]*mfromkm^2;
% Mean-squared differentiability
nu=[2 2 3 3];
% Correlation parameter [m]
% This needs to alternate so the x-axes are common
rho=[50 100 50 100]*mfromkm;

% Number of simulations per processor, higher is better for small grids
npr=5;

% Y axis limits
yls=[-0.1 1.1];

% Sizes of the fields under investigation
% Must be able to quarter these!
Ns=8:2:128;
% Physical unit dimensions
p.dydx=[10 10]*mfromkm;
% Quartering OFF for EGGERS1
p.quart=0;
% Up the blurring! Odds are always good since they turn evens and into
% evens and keep odds being odds. This works better than before,
% especially for the low numbers, very noticeably! Hashes saved.
% p.blurs=3;
% Actually, use BLUROSY. Hashes saved.
p.blurs=-1;

% Four panels with different parameters
clf
[ah,ha,H]=krijetem(subnum(2,2));

% For each of the correlation lengths
for fndex=1:length(rho)
  % Supply the Matern parameters
  th0=[s2 nu(fndex) rho(fndex)];
  % Make this a hash, too!

  % For this set of parameters, make a unique hashed filename
  fname=hash([struct2array(orderfields(p)) th0 Ns npr matlabpool('size')],'SHA-1');
  fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('EGGERS1_%s.mat',fname));

  if ~exist(fnams,'file') 
    % For each of the data sizes
    for lndex=1:length(Ns)
      p.NyNx=[Ns(lndex) Ns(lndex)];
      % Divide the workload over the processors; note: you must matlabpool!
      spmd
        for index=1:npr
          % Perform a simulation
          [Hx,th0p,pp]=simulosl(th0,p); 
          % Naive estimation of the variance
          vHx(index)=var(Hx);
        end
      end
      % Combine the results from all the processors
      nvars=cat(2,vHx{:});
      % The mean of these guys
      mvars(lndex)=mean(nvars);
      % Collapse the parameters, keep track of physical lengths
      pp=pp{1}; 
      plen(lndex)=sqrt(prod(pp.NyNx))*sqrt(prod(pp.dydx));
      % Bias of the naive estimator using observed covariance 
      b1(lndex)=varbias(th0,pp,1,0);
      % Bias of the naive estimator using blurred likelihood 
      b3(lndex)=varbias(th0,pp,3,0);
      % Bias of the naive estimator using analytical full likelihood 
      b4(lndex)=varbias(th0,pp,4,0);
    end
    save(fnams,'plen','mvars','th0','b1','b3','b4','yls','mfromkm','p','fndex')
  else
    disp(sprintf('Loading %s',fnams))
    load(fnams)
  end
        
  % Make the plot
  axes(ah(fndex))
  % Order from 'worst' to 'best' for display  
  hold on

  % Plot the correlation length times pi times 2
  plc(fndex)=plot(2*pi*[th0(3) th0(3)]/mfromkm,yls,'k-');
  % Prediction using the analytical full likelihood
  pl(fndex,4)=plot(plen/mfromkm,[th0(1)-b4]/th0(1),'r-'); 
  % Prediction using the blurred likelihood
  pl(fndex,3)=plot(plen/mfromkm,[th0(1)-b3]/th0(1),'k-');
  % Prediction using the full summation of the covariance
  pl(fndex,2)=plot(plen/mfromkm,[th0(1)-b1]/th0(1),'b-');
  % Mean of the poor estimates, as actually simulated and averaged
  pl(fndex,1)=plot(plen/mfromkm,mvars/th0(1),     'ko'); 
    
  hold off

  % Cosmetics
  ylim(yls)
  xl(fndex)=xlabel(sprintf('grid size (km) ; %s = %i pixels',...
                           '\pi\rho',round(pi*th0(3)/sqrt(prod(p.dydx)))));
  % yl(fndex)=ylabel(sprintf('normalized expected sample variance'));
  yl(fndex)=ylabel(sprintf('observed to predicted <s^2>/%s^2','\sigma'));
  tl(fndex)=title(sprintf('%s = %3.1f km ; %s = %i ; %s = %i km',...
                   '\sigma',sqrt(th0(1))/mfromkm,...
                   '\nu',th0(2),...
                   '\rho',th0(3)/mfromkm));
  
  set(pl(:,2:end),'linew',1)
  
  if fndex==1
    leg=legend(pl,...
               sprintf('mean (%i%s)',...
                       npr*matlabpool('size'),'\times'),...
               'full-covariance',...
               'blurred-likelihood',...
               'full-likelihood',...
               'Location','SouthEast');
  end
  drawnow 
end

% Cosmetics
longticks(ah)
set(pl(:,1),'Color','k','Marker','o','MarkerFaceC','w','MarkerEdgeC','k','MarkerS',4)
set(ah,'ygrid','on','box','on')
set(tl,'FontS',12)
nolabels(ha([3 4]),2); delete(yl([2 4]))
nolabels(ah([1 2]),1); delete(xl([1 2]))
serre(H,[],'across'); serre(H',1/3,'down')
legend

% Printo
fig2print(gcf,'portrait')
figna=figdisp([],[],[],1);
system(sprintf('epstopdf %s.eps',figna));
system(sprintf('rm -f %s.eps',figna));

