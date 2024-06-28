function fit = BIClassoVAR(Ydata,q,intr,normalize,tol_glmnet)

if nargin < 2 || isempty(q)
    q=1;        % default to VAR(1) fit
end
if nargin < 3 || isempty(intr)
    intr=1;     % default is to include intercepts (unpenalized)
end
if nargin < 4 || isempty(normalize)
    normalize=0;    % default is not to normalize regressors
end
if nargin < 5 || isempty(tol_glmnet)
    tol_glmnet=1e-4;   % glmnet tolerance
end

% Unpack data
[nplusq,p]=size(Ydata);     % data dimensions
n=nplusq-q;                 % effective sample size with q lags
qp=q*p;                     % #parameters
Y=Ydata(q+1:q+n,:);         % outcomes
X=nan(n,qp);                % regressors
for ell=1:q
    block_ell=(ell-1)*p+1:ell*p;     % where to store ellth lag
    lag_ell=q+1-ell:q+n-ell;         % periods corresponding to ellth lag
    X(:,block_ell)=Ydata(lag_ell,:); % store 1st lags first, 2nd second,...
end

if intr
    opts.intr=true;             % fit an intercept
    const=nan(p,1);             % intercept placeholder
else
    opts.intr=false;            % o/w don't
    const=[];
end
if normalize
    opts.standardize=true;      % do standardize (if on different scales)
else
    opts.standardize=false;     % don't standardize (presumed equal scale)
end
opts.thresh=tol_glmnet;         % tolerance for coordinate descent
options=glmnetSet(opts);        % apply options for glmnet

That=nan(p,qp);                 % slopes placeholder
for i=1:p
    bic_i=BICglmnet(X,Y(:,i),options);
    That(i,:)=bic_i.beta;
    if intr
        const(i)=bic_i.a0;
    end
end

fit.const=const;                % intercepts (if any)
fit.That=That;                  % slopes

end