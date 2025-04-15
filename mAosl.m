function varargout=mAosl(k,th,params,xver)
% [mth,mththp,A,k]=mAosl(k,th,params,xver)
%
% Calculates auxiliary derivative quantities for Olhede & Simons (2013)
% in the SINGLE-FIELD case. With some obvious redundancy for A.
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with elements [not scaled]:
%          th(1)=s2  The first Matern parameter [variance in unit^2]
%          th(2)=nu  The second Matern parameter [differentiability]
%          th(3)=rh  The third Matern parameter [range in m]
% params   blurs 0 No wavenumber blurring
%                1 No wavenumber blurring, effectively
%                N Convolutional blurring, errors 
%               -1 Exact blurring
%              Inf Exact blurring, effectively
% xver     Excessive verification [0 or 1]
%
% OUTPUT:
%
% mth      The first-partials "means" for GAMMIOSL, HESSIOSL, as a cell
% mththp   The second-partials parameters for GAMMIOSL, HESSIOSL, as a cell
% A        The "A" matrices for GAMMIOSL, HESSIOSL, as a cell
% k        The actual wavenumbers being returned (could be different than requested
%          if blurring is being performed, which recomputed k)
%
% Last modified by fjsimons-at-alum.mit.edu, 04/15/2025
% Last modified by olwalbert-at-princeton.edu, 04/15/2025

% Extra verification?
defval('xver',1)

% The number of parameters to solve for
np=length(th);

% If we didn't ask for blurring, make sure it's explicitly unblurred
defstruct('params',{'blurs'},0)

% Extract the parameters from the input
s2=th(1);
nu=th(2);
rh=th(3);

% If somehow it was called with Inf we meant exact blurring
% in the sense of MATERNOSP and BLUROSY
if isinf(params.blurs); params.blurs=-1; end

% We are able to blur, exactly
if params.blurs==-1
    % We calculate the derivatives of the blurred spectral density from the 
    % blurred autocovariance via MATERNOSY which is called by SPATMAT
    [dSbards2,kk]=blurosy(th,params,[],[],[],1);
    dSbardnu=blurosy(th,params,[],[],[],2);
    dSbardrh=blurosy(th,params,[],[],[],3);
    % Get the exactly blurred spectrum
    Sbar=maternosp(th,params);
    % Use eq. (112) in doi: 10.1093/gji/ggt056
    mth{1}=dSbards2./Sbar;
    mth{2}=dSbardnu./Sbar;
    mth{3}=dSbardrh./Sbar;

    % If we'd made too many...
    if prod(size(k))~=prod(size(kk)) && all(kk(~~kk)==k)
        % ... cut them to the same size
        for ind=1:np
            mth{ind}=mth{ind}(~~kk);
        end
    end
    
    if nargout>1
        % We will have to build in second partials for mththp and A, e.g. for HESSIOSL
        % but we won't be needing them, and they are also never needed for FISHIOSL
        % warning('Returning UNBLURRED values for calculations of mththp and A')
        params.blurs=0;
        [~,mththp,A]=mAosl(k,th,params,xver);
    else
        xver=0;
        mththp=NaN;
        A=NaN;
    end
else
    if isinf(nu)
        % When we are considering the special case of the limit of the spectral
        % density as nu approaches infinity, we must calculate the first and second
        % derivatives from the limit of the spectral density
        % ms2
        mth{1}=1/s2;
        % mnu 
        mth{2}=0;
        % for d-dimensions, mth{3} = d/rh - pi^2*rh*k(:).^2/2
        % mrh for 2-dimensions
        mth{3}=2/rh-pi^2*rh*k(:).^2/2;
        
        if nargout>1
            % dms2/ds2
            mththp{1}=-1/s2^2;
            % dmnu/dnu
            mththp{2}=0;
            % for d-dimensions, mththp{3} = -d/rh^2-pi^2*k(:).^2/2
            % dmrh/drh for 2-dimensions
            mththp{3}=-2/rh^2-pi^2*k(:).^2/2;
            % Cross-terms
            mththp{4}=0;
            mththp{5}=0;
            % dmrhdnu or dmnudrh
            mththp{6}=0;
        else
            mththp=NaN;
        end
    else
        % Auxiliary variable
        avark=4*nu/pi^2/rh^2+k(:).^2;

        % The "means", which still depend on the wavenumbers, sometimes
        % The first derivatives of the logarithmic spectral density
        % Eq. (A25) in doi: 10.1093/gji/ggt056
        mth{1}=1/s2;
        % Eq. (A26) in doi: 10.1093/gji/ggt056
        mth{2}=(nu+1)/nu+log(4*nu/pi^2/rh^2)...
               -4*(nu+1)/pi^2/rh^2./avark-log(avark);
        % Eq. (A27) in doi: 10.1093/gji/ggt056
        mth{3}=-2*nu/rh+8*nu/rh*(nu+1)/pi^2/rh^2./avark;

        if nargout>1
            % The second derivatives of the logarithmic spectral density
            vpiro=4/pi^2/rh^2;
            avark=nu*vpiro+k(:).^2;
            % Here is the matrix of derivatives of m, diagonals first, checked with 
            % syms s2 nu rh k avark vpiro pi
            % Checked dms2/ds2
            mththp{1}=-1/s2^2;
            % Checked dmnu/dnu
            mththp{2}=2/nu-(nu+1)/nu^2-2*vpiro./avark+(nu+1)*vpiro^2./avark.^2;
            % Checked dmrh/drh
            mththp{3}=2*nu*(1-3*(nu+1)*vpiro./avark+2*nu*(nu+1)*vpiro^2./avark.^2)/rh^2;
            % Here the cross-terms, which are verifiably symmetric
            mththp{4}=0;
            mththp{5}=0;
            % This I've done twice to check the symmetry, checked dmrhdnu or dmnudrh
            mththp{6}=2/rh*(-1+[2*nu+1]*vpiro./avark-nu*[nu+1]*vpiro^2./avark.^2);
        else 
            mththp=NaN;
        end
    end
    kk=k;
end

% Full output and extra verification
if nargout>2 || xver==1
  % How many wavenumbers?
  lk=length(kk(:));
  
  % Initialize
  A=cellnan(length(th),lk,3);

  % Eq. (A28) in doi: 10.1093/gji/ggt056
  A{1}=-repmat(mth{1},lk,3);
  % Eq. (A29) in doi: 10.1093/gji/ggt056
  if nu==Inf
    A{2}=-repmat(mth{2},lk,3);
  else
      A{2}=-repmat(mth{2},1,3);
  end
  % Eq. (A30) in doi: 10.1093/gji/ggt056
  A{3}=-repmat(mth{3},1,3);

  % Verification mode
  if xver==1
    % Inf case used to be excluded, until TRACECHECK was revisited
    % In this univariate case, we have some rather simple forms
    % In the bivariate case, these are the nonzero lower-triangular
    % entries of a matrix that is another's Cholesky decomposition
    L=[ones(lk,1) zeros(lk,1) ones(lk,1)];
    tracecheck(L,A,mth,9,1)
    % How about the second TRACECHECK? We should be able to involve mththp
  end
else
  A=NaN;
end

% As much output as you desire
varns={mth,mththp,A,kk};
varargout=varns(1:nargout);
