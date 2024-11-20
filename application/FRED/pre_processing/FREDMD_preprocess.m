% Initiate
clc;clear;format compact;

% Path
cd 'C:\Users\kzb125\Downloads\hd-var\application\FRED\pre_processing'

% -------------------------------------------------------------------------
% File name
csv_in='FRED-MD_2022-05.csv';

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

out_folder = 'C:\Users\kzb125\Downloads\hd-var\application\FRED\data\';
writematrix(data_imp,sprintf('%sFRED-MD_2022-05_preprocessed.csv', out_folder))
