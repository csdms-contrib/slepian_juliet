function eggers1
% EGGERS1
%
% Makes a FIGURE illustrating, by simulation (using SIMULOSL), and
% predicting, via three different analytical approximations (using VARBIAS),
% the biasing influence that unknown spatial covariance (primarily via the
% differentiability and range parameters nu and rho) has on the brute-force
% estimation of the variance (sigma^2) of an isotropic Matern process.
%
% SEE ALSO:
%
% SIMULOSL, VARBIAS
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu 06/23/2018

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

% Experimental parameters, in standard units
% Sizes of the fields under investigation (watch if p.quart is on!)
Ns=8:2:128;
% Physical unit dimensions
p.dydx=[10 10]*mfromkm;
% Quartering (off, p.quart=0, for the published figure)
p.quart=0;
% Blurring (using BLUROSY, p.blurs=-1, for the published figure)
% If using the inexact regridded convolutional BLUROS procedure,
% relatively large, positive, odd numbers are a good choice, since they
% keep the parity of even- and odd-dimensional grids 
p.blurs=-1;

% Number of simulations per processor, higher is better for small grids
npr=5;
% Number of processors
NumWorkers=8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the figure: four panels with different parameters
clf; [ah,ha,H]=krijetem(subnum(2,2));

% Y axis limits
yls=[-0.1 1.1];

% For each of the correlation lengths
for fndex=1:length(rho)
  % Supply the Matern parameters
  th0=[s2 nu(fndex) rho(fndex)];

  % For this ordered (!) set of parameters, make a unique hashed filename
  % so that you can save the results and don't have to redo the computations
  fname=hash([struct2array(orderfields(p)) th0 Ns npr NumWorkers],'SHA-1');
  % You need to have an environmental variable file structure set up
  fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s_%s.mat',upper(mfilename),fname));
  
  % If it hasn't been precomputed and saved ahead of time
  if ~exist(fnams,'file')
    % Initialize the pool of workers
    if isempty(gcp('nocreate')); pnw=parpool(NumWorkers); end
    % For each of the data sizes
    for lndex=1:length(Ns)
      p.NyNx=[Ns(lndex) Ns(lndex)];
      % Divide the workload over the processors
      spmd
        for index=1:npr
          % Perform a simulation
          [Hx,th0p,pp]=simulosl(th0,p); 
          % Naive estimation of variance: this is what we are debunking
          vHx(index)=var(Hx);
        end
      end
      % Combine the results from all the processors
      nvars=cat(2,vHx{:});
      % The mean of all the naive variance estimates
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
    % Close the pool of workers
    delete(pnw)
    % Save into the hash so the above won't need to be recalculated next time
    save(fnams,'plen','mvars','th0','b1','b3','b4','yls','mfromkm','p','fndex')
  else
    disp(sprintf('%s loading %s',upper(mfilename),fnams))
    load(fnams)
  end
        
  % Make the plot
  axes(ah(fndex))
  % Order from 'worst' to 'best' for display  
  hold on

  % Plot empirical distance where bias decreases to about a third of the truth
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
  % yl(fndex)=ylabel(sprintf('observed to predicted <s^2>/%s^2','\sigma'));
  yl(fndex)=ylabel(sprintf('%s estimates relative to truth','\sigma^2'));
  tl(fndex)=title(sprintf('%s = %3.1f km ; %s = %i ; %s = %i km',...
                   '\sigma',sqrt(th0(1))/mfromkm,...
                   '\nu',th0(2),...
                   '\rho',th0(3)/mfromkm));
  
  set(pl(:,2:end),'linew',1)
  
  if fndex==1
    leg=legend(pl,...
               sprintf('mean (%i%s)',...
                       npr*NumWorkers,'\times'),...
               'full-covariance',...
               'blurred-likelihood',...
               'full-likelihood',...
               'Location','SouthEast');
  end
  drawnow 
end

% Cosmetics to match EGGERS4 exactly
longticks(ah)
set(pl(:,1),'Color','k','Marker','o','MarkerFaceC','w','MarkerEdgeC','k','MarkerS',4)
set(ah,'ygrid','on','box','on')
set(tl,'FontSize',12)
nolabels(ha([3 4]),2); delete(yl([2 4]))
nolabels(ah([1 2]),1); delete(xl([1 2]))
serre(H,[],'across'); serre(H',1/3,'down')
legend;

% Prin to file
fig2print(gcf,'portrait')
figdisp([],[],[],2,'epsc','epstopdf');
