function b=varbias(th0,p,calx,calxo)
% b=VARBIAS(th0,p,calx,calxo)
%
% Computes the (NEGATIVE) bias affecting the SAMPLE variance of a spatially
% CORRELATED random Matern field when it is estimated while ignoring this
% correlation (because it too, needs to be estimated!) So... if you simulate
% some data using Hx=SIMULOSL(p,th0), and you simply compute s2hat=var(Hx),
% the expectation of this estimate is the true variance, s2, indeed th0(1),
% MINUS the bias computed here using a variety of ways.
%
% INPUT:
%
% th0      The three parameters of the (Matern) covariance model
% p        A parameter vector with at least:
%          NyNx size of the data set
%          dydx physical units of the data set
% calx     1 full sum of all spatial covariance pairs [default]
%          2 full sum of all spatial covariance pairs, slower method
%          3 approximation using the blurred likelihood, very good
%          4 approximation using the analytic full likelihood; this is
%             lightningly fast but only applies when the integral of the
%             correlation function can be well approximated by the Riemann
%             sum, which is rare 
% calxo    1 override selection to speedier calx=3 [default, see EGGERS4]
%          0 no override [see EGGERS1, which illustrates the methods]
% 
% OUTPUT:
%
% b        The bias! 
%
% SEE ALSO:
%
% SIMULOSL
%
% EXAMPLE:
%
% for ind=1:40
%   p.NyNx=[64 64];
%   [Hx,th0,p]=simulosl([],p); vHx(ind)=var(Hx);
% end
% b=varbias(th0,p);
% [mean(vHx) th0(1)-b]
%
% EXAMPLE: See usage in EGGERS1, EGGERS3, and so on 
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 12/02/2016

defval('xver',0)

% Sadly, essentially, you really NEVER want to do option 2
defval('calx',1)
defval('calxo',1)
% Default Matern parameters
defval('th0',[100000 2 50000])
% Supply the needed parameters, keep the givens, extract to variables
fields={               'dydx','NyNx','blurs','kiso','quart'};
defstruct('p',fields,...
	    {                      [20 20]*1e3,[64 64],2,NaN,1});
% Flatten the structure for convenience
struct2var(p)

% If the fields are BIG enough you can just switch to calx=3!  See
% EGGERS4... but obviously cannot do this for EGGERS1, as this is testing
% the very approximation under which this may hold or not
if calxo==1
  if sqrt(prod(NyNx))*sqrt(prod(dydx))>4*th0(3);
    calx=3;
    disp(sprintf('%s field big enough to switch to calx = %i',upper(mfilename),calx));
  end
end
% In case you've done any of this before... ORDER matters
fname=hash([struct2array(orderfields(p)) th0 calx],'SHA-1');
fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s_%s.mat',upper(mfilename),fname));

if ~exist(fnams,'file')
  t=tic;
  if calx==1
    % All the distances from anywhere to anywhere
    Dxxp=xxpdist([0:NyNx(1)-1]*dydx(2),[0:NyNx(2)-1]*dydx(2));
    
    % Convert to covariance
    Cy=maternosy(Dxxp,th0);
    % Sum and add all at the same time
    b=sum(sum(Cy))/prod(NyNx)^2;
  elseif calx==2
    % First, all the distance to the first grid point
    D1xp=sspdist(1,1:prod(NyNx),NyNx,dydx)';
    % Maybe extract the unique ones from it to save time?
    % Convert to covariance
    Cy=maternosy(D1xp,th0);
    % Initialize
    b=0;
    for sndex=2:prod(NyNx)
      % Compute the other times these values are being hit
      % Maybe adjust the loop inside to only do the right half?
      [~,pm]=ssp21sp(sndex,NyNx);
      % At any rate, you need to sum and add
      b=b+sum(Cy(pm(1,:)).*pm(2,:));
    end
    b=b/prod(NyNx)^2;
  elseif calx==3
    % Find the zero wavenumber
    kzero=sub2ind(NyNx,floor(NyNx(1)/2)+1,floor(NyNx(2)/2)+1);
    % Let's use the blurred spectral density
    Sb=maternosp(th0,p,xver);
    % Pick out the zero-wavenumber blurred version... should be very good
    b=Sb(kzero)/prod(NyNx)/prod(dydx)*(2*pi)^2;
  elseif calx==4
    % The first line is identical to the second line below
    b1=maternos(0,th0)/prod(NyNx)/prod(dydx)*(2*pi)^2;
    b=th0(1)/prod(NyNx)/prod(dydx)*pi^3*th0(3)^2;
    % Check for consistency but don't make a big deal out of it
    difer(b1-b,8,[],NaN)
    % But either of them relies on dense sampling relative to rho so this
    % won't be always applicable
  end
  disp(sprintf('%s took %f seconds',upper(mfilename),toc(t)))
  save(fnams,'b','p','th0','calx')
else
  disp(sprintf('%s loading %s',upper(mfilename),fnams))
  load(fnams)
end
