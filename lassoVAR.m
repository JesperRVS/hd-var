function fit = lassoVAR(Ydata,q,post,intr,c,gamma,K,tol_Ups,tol_glmnet,nowarn)

if nargin < 2 || isempty(q)
    q=1;        % default to VAR(1) fit
end
if nargin < 3 || isempty(post)
    post=1;     % default to refitting post selection (in every step)
end
if nargin < 4 || isempty(intr)
    intr=1;     % default is to include intercepts (unpenalized)
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

% Penalty level
if nargin < 5 || isempty(c)
   c=1.1; % default markup is 10 pct
end
mstar=floor(n^(1/5));   % deviation bound optimal block size
lstar=floor(n/mstar);   % implied #blocks and...
if nargin < 6 || isempty(gamma)
    gamma=.1/log(max(n,qp));   % default probability tolerance
end
% penalty implied by deviation bound optimal block size:
lambdastar=((2*c*n)/sqrt(mstar*lstar))*norminv(1-gamma/(2*q*p^2));
lambda_glmnet=lambdastar/(2*n); % penalty in eyes of glmnet
% note: glmnet defines "lambda" based on *half* of *average* square loss,
% while our "lambda" stems from *sum* square loss (w/o the 1/2)

% Set tolerances
if nargin < 7 || isempty(K)
    K=15;       % at most K updates
end
if nargin < 8 || isempty(tol_Ups)
    tol_Ups=1e-3;   % loading update tolerance (relative change so .1%)
end
if nargin < 9 || isempty(tol_glmnet)
    tol_glmnet=1e-4;    % glmnet tolerance (their default)
end
if nargin < 10 || isempty(nowarn)
    nowarn=0;   % default is to report lack of update convergence, if any
else
    nowarn=1;
end

if intr % if intercepts are included, demean before continuing
    Ybar=mean(Y,1);
    Y=Y-Ybar;
    Xbar=mean(X,1);
    X=X-Xbar;
end
% Note: Means stored to back out intercept estimates

% Initial step
UpsInit=sqrt((Y.^2)'*(X.^2)/n);                         % loadings
ThatInit=mlasso(X,Y,lambda_glmnet,UpsInit,tol_glmnet);  % estimates
if post
    sel=(ThatInit~=0);                              % initial selection
    for i=1:p 
        ThatInit(i,sel(i,:))=X(:,sel(i,:))\Y(:,i);  % refit post selection
    end
end
if intr                             % if intercepts requested...
    constInit=Ybar'-ThatInit*Xbar'; % ...back them out
else
    constInit=[];
end

% At most K loading updates using LASSO in every step
UpsRefi=nan(p,qp,K);
ThatRefi=nan(p,qp,K);
if intr                 % if intercepts requested...
    constRefi=nan(p,K); % ...create placeholder
else
    constRefi=[];
end
for k=1:K
    if k==1 % using initial estimates to calculate residuals
        Res_kminus1=Y-X*ThatInit';
    else    % using updated estimates to calculate residuals
        Res_kminus1=Y-X*ThatRefi(:,:,k-1)';
    end
    UpsRefi(:,:,k)=sqrt((Res_kminus1.^2)'*(X.^2)/n);
    ThatRefi_k=mlasso(X,Y,lambda_glmnet,UpsRefi(:,:,k),tol_glmnet);
    if post % refitting post selection
        sel=(ThatRefi_k~=0);                                % selection
        for i=1:p 
            ThatRefi_k(i,sel(i,:))=X(:,sel(i,:))\Y(:,i);    % refit
        end
    end
    ThatRefi(:,:,k)=ThatRefi_k;
    if intr                                    % if intercepts requested...
        constRefi(:,k)=Ybar'-ThatRefi_k*Xbar'; % ...back them out
    end
    % Check if penalty loading tolerance is met
    if k==1
        UpsOld=UpsInit;
    else
        UpsOld=UpsRefi(:,:,k-1);
    end
    UpsNew=UpsRefi(:,:,k);
    dUps=UpsNew-UpsOld;
    reldiffUps=norm(dUps(:))/norm(UpsOld(:));   % if change in loadings...
    if reldiffUps<=tol_Ups                      % ...falls below tolerance...
        break;                                  % ...stop.
    end                                         % o/w continue updating
end
Ups=UpsRefi(:,:,k);         % loadings after updating
That=ThatRefi(:,:,k);       % slope estimates after updating
if intr                     % if intercepts requested...
    const=Ybar'-That*Xbar'; % ...back them out.
else
    const=[];               % otherwise leave empty
end

if nowarn==0
    if k==K && reldiffUps>tol_Ups % warn if tolerance not met within K updates
    warning('Maximum number of updates reached. Consider increasing K.')
    fprintf(sprintf('Relative change in penalty loadings is %3.1g percent\n',100*reldiffUps))
    end
end

fit.lambda=lambdastar;      % penalty level
fit.UpsInit=UpsInit;        % initial penalty loadings
fit.UpsRefi=UpsRefi;        % sequence of refined penalty loadings
fit.Ups=Ups;                % final penalty loadings
fit.constInit=constInit;    % initial intercept
fit.constRefi=constRefi;    % sequence of intercepts from refinements
fit.const=const;            % final intercept
fit.ThatInit=ThatInit;      % initial slopes
fit.ThatRefi=ThatRefi;      % sequence of slopes from refinements
fit.That=That;              % final slopes
fit.kterm=k;                % #updates used
fit.reldiffUps=reldiffUps;  % final relative change in penalty loadings

end