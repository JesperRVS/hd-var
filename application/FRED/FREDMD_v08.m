% Initiate
clc;clear;format compact;

% Path
cd 'C:\Users\kzb125\Downloads\hd-var\application\FRED'

% -------------------------------------------------------------------------
% File name
csv_in=readmatrix('data/FRED-MD_2022-05.csv');

% =========================================================================
% PART 1: LOAD AND LABEL DATA (from McCracken & Ng)
dum=importdata(csv_in,',');     % load data
series=dum.textdata(1,2:end);   % variable names
tcode=dum.data(1,:);            % transformation numbers
rawdata=dum.data(2:end,:);      % raw data

% Month/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_month=final_datevec(2);
final_year=final_datevec(1);

% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
dates = (1959+1/12:1/12:final_year+final_month/12)';

% T = number of months in sample
T=size(dates,1);
rawdata=rawdata(1:T,:);

% =========================================================================
% PART 2: PROCESS DATA (from McCracken & Ng)

% Transform raw data to be stationary using auxiliary function
% prepare_missing()
yt=prepare_missing(rawdata,tcode);

% Reduce sample to usable dates: remove first two months because some
% series have been first differenced
yt=yt(3:T,:);
dates=dates(3:T,:);

% Remove outliers using auxiliary function remove_outliers(); see function
% or readme.txt for definition of outliers
%   data = matrix of transformed series with outliers removed
%   n = number of outliers removed from each series
[data,~]=remove_outliers(yt);

% =========================================================================
% PART 3: IMPUTATION (from initialization of McCracken & Ng's EM algorithm)
% Fill in missing values for each series with the average of that series.
data_imp=data;  % copy and...
for i=1:numel(series) 
    data_imp(isnan(data(:,i)),i)=mean(data(:,i),'omitnan'); % ...impute
end
% Note: noise added due to imputation abstracted from in our theory.

% =========================================================================
% PART 4: 1-PERIOD-AHEAD FORECASTING VIA VAR(q) LASSO AND POST-LASSO
[n_all,p]=size(data_imp);   % #time periods total
% qlist=1:12;                 % lags
qlist=1:4;                 
Q=numel(qlist);             % #lags
% N=120;                      % #1-period-ahead forecasts to be made (10 yrs)
N=12;                      
qmax=max(qlist);            % maximal lag
n=n_all-N-qmax;             % implied effective sample size for forecasting

Lassos    =cell(N,Q);
Posts     =cell(N,Q);
BICs_raw  =cell(N,Q);
BICs_norm =cell(N,Q);
ARs_consts=cell(N,Q);
ARs_slopes=cell(N,Q);
OLS_consts=cell(N,Q);
OLS_slopes=cell(N,Q);
OLS_illcond=nan(N,Q);
tic;
for thisq=1:Q
    q=qlist(thisq);
    fprintf(sprintf('Forecasting initiated for q=%d\n',q));
%     for thisnplus1=1:N
    parfor thisnplus1=1:N
        nplus1=qmax+n+thisnplus1;
        sample=nplus1-(n+q):nplus1-1;   % effective sample of size n
        YEst=data_imp(sample,:);        % ... prior to n+1
        % VAR(q) Lasso (w/ const)
        Lassos{thisnplus1,thisq}=lassoVAR(YEst,q,0,1,[],[],[],[],[],1); 
        % VAR(q) Post-Lasso (w/ const)
        Posts{thisnplus1,thisq} =lassoVAR(YEst,q,1,1,[],[],[],[],[],1); 
        % VAR(q) BIC-Lasso w/ raw data (w/ const)
        BICs_raw{thisnplus1,thisq}=BIClassoVAR(YEst,q,1,0);
        % VAR(q) BIC-Lasso w/ standardization (w/ const)
        BICs_std{thisnplus1,thisq}=BIClassoVAR(YEst,q,1,1);
        % AR(q) OLS (w/ const)
        Y=YEst(q+1:q+n,:);              % outcomes
        qp=q*p;                         % #parameters
        const_temp=nan(p,1);            % intercepts p x 1
        Theta_temp=nan(p,qp);           % slopes, storing as p x qp matrix
        for i=1:p
            X_i=nan(n,q);
            for ell=1:q
                X_i(:,ell)=YEst(q+1-ell:q+n-ell,i);
            end
            LS_i=[ones(n,1),X_i]\Y(:,i);
            const_temp(i)=LS_i(1);
            AR_entries_i=i:p:qp;
            Theta_temp(i,AR_entries_i)=LS_i(2:1+q); 
            nonAR_entries_i=setdiff(1:qp,i:p:qp);
            Theta_temp(i,nonAR_entries_i)=0;
        end
        ARs_consts{thisnplus1,thisq}=const_temp;
        ARs_slopes{thisnplus1,thisq}=Theta_temp;
        % VAR(q) OLS (w/ const)
        X=nan(n,qp);                    % regressors
        for ell=1:q
            block_ell=(ell-1)*p+1:ell*p;     % where to store ellth lag
            lag_ell=q+1-ell:q+n-ell;         % periods corresponding to ellth lag
            X(:,block_ell)=YEst(lag_ell,:);  % store 1st lags first, 2nd second,...
        end
        OLS_illcond(thisnplus1,thisq)=isIllConditioned(decomposition(X));
        if ~OLS_illcond(thisnplus1,thisq)
            LS=[ones(n,1),X]\Y;
            OLS_consts{thisnplus1,thisq}=LS(1,:)';
            OLS_slopes{thisnplus1,thisq}=LS(2:end,:)';
        end
    end
    fprintf(sprintf('Done with q=%d\n',q));
    toc;
end
fprintf('Forecasting completed.\n');
toc;

InvSampleVar=diag(1./var(data_imp));   % inverse sample variances (for weighting)
WeightedSFEl2_lasso=nan(N,Q);          % placeholders, weighted ell_2 FEs
WeightedSFEl2_post =nan(N,Q);
WeightedSFEl2_bic_raw  =nan(N,Q);
WeightedSFEl2_bic_std  =nan(N,Q);
WeightedSFEl2_ar   =nan(N,Q);
WeightedSFEl2_ols  =nan(N,Q);
for thisq=1:Q
    q=qlist(thisq);
    for thisnplus1=1:N
        nplus1=qmax+n+thisnplus1;
        Ynplus1=data_imp(nplus1,:)';    % Y_{n+1}, to be forecast
        YPre=nan(q,p);                  % predictors Y_{n},...,Y{n-q}
        for ell=1:q                     % ...as q x p matrix
            YPre(ell,:)=data_imp(nplus1-ell,:); % going back in time
        end
        YPre=YPre';         % ... as p x q matrix
        YPre=YPre(:);       % ... as pq x 1 vector
        % VAR(q) Lasso
        Ynplus1hat_lasso=Lassos{thisnplus1,thisq}.const...
            +Lassos{thisnplus1,thisq}.That*YPre;  % 1-period-ahead forecast
        FE_lasso=Ynplus1hat_lasso-Ynplus1;                  % ...error vector
        WeightedSFEl2_lasso(thisnplus1,thisq)=FE_lasso'*InvSampleVar*FE_lasso;
        % VAR(q) Post-Lasso
        Ynplus1hat_post=Posts{thisnplus1,thisq}.const....
            +Posts{thisnplus1,thisq}.That*YPre;  % 1-period-ahead forecast
        FE_post=Ynplus1hat_post-Ynplus1;                  % ...error vector
        WeightedSFEl2_post(thisnplus1,thisq)=FE_post'*InvSampleVar*FE_post;
        % VAR(q) BIC-Lasso w/ raw data
        Ynplus1hat_bic_raw=BICs_raw{thisnplus1,thisq}.const....
            +BICs_raw{thisnplus1,thisq}.That*YPre;  % 1-period-ahead forecast
        FE_bic_raw=Ynplus1hat_bic_raw-Ynplus1;                  % ...error vector
        WeightedSFEl2_bic_raw(thisnplus1,thisq)=FE_bic_raw'*InvSampleVar*FE_bic_raw;
        % VAR(q) BIC-Lasso w/ standardization
        Ynplus1hat_bic_std=BICs_std{thisnplus1,thisq}.const....
            +BICs_std{thisnplus1,thisq}.That*YPre;  % 1-period-ahead forecast
        FE_bic_std=Ynplus1hat_bic_std-Ynplus1;                  % ...error vector
        WeightedSFEl2_bic_std(thisnplus1,thisq)=FE_bic_std'*InvSampleVar*FE_bic_std;
        % AR(q) OLS
        Ynplus1hat_ar=ARs_consts{thisnplus1,thisq}....
            +ARs_slopes{thisnplus1,thisq}*YPre;  % 1-period-ahead forecast
        FE_ar=Ynplus1hat_ar-Ynplus1;                  % ...error vector
        WeightedSFEl2_ar(thisnplus1,thisq)=FE_ar'*InvSampleVar*FE_ar;
        if ~OLS_illcond(thisnplus1,thisq)
            % VAR(1) OLS
            Ynplus1hat_ols=OLS_consts{thisnplus1,thisq}...
                +OLS_slopes{thisnplus1,thisq}*YPre;  % 1-period-ahead forecast
            FE_ols=Ynplus1hat_ols-Ynplus1;                  % ...error vector
            WeightedSFEl2_ols(thisnplus1,thisq)=FE_ols'*InvSampleVar*FE_ols;
        end
    end
end

save(sprintf('%s_N%d_Q%d',date,N,Q),...
     'WeightedSFEl2_lasso','WeightedSFEl2_post','WeightedSFEl2_bic_raw',...
     'WeightedSFEl2_bic_std','WeightedSFEl2_ar','WeightedSFEl2_ols');

% save(sprintf('ApplicationResults_%s_N%d_Q%d', date, N, Q));
if isunix && ~ismac
    quit;
end


%%
