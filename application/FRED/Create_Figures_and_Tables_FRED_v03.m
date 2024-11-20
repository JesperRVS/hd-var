% Initiate
clc;clear;format compact;

% Set directory and paths
if ispc
    % Paths when working from office
    cd 'D:\Dropbox\High-Dimensional Time Series\Matlab\hdvar\application\FRED';
    addpath 'D:\Dropbox\High-Dimensional Time Series\Matlab\hdvar';
    addpath 'D:\Dropbox\High-Dimensional Time Series\Matlab';
elseif isunix
    if ismac
        % Paths when working from home
        cd '/Users/jesperrvs/Library/CloudStorage/Dropbox/High-Dimensional Time Series/Matlab/hdvar/application/FRED';
        addpath '/Users/jesperrvs/Library/CloudStorage/Dropbox/High-Dimensional Time Series/Matlab/hdvar';
        addpath '/Users/jesperrvs/Library/CloudStorage/Dropbox/High-Dimensional Time Series/Matlab';
        addpath '/Users/jesperrvs/Library/CloudStorage/Dropbox/High-Dimensional Time Series/Matlab/glmnet-matlab-master';
    end
end

%% Load results
% load('2024Feb08_N120_Q12','-mat');
load('08-Feb-2024_N120_Q12','-mat');
% load('09-Feb-2024_N12_Q4','-mat');

%%
figure;
VAR1_lasso=mean(WeightedSFEl2_lasso(:,1));
plot(1:12,100*[mean(WeightedSFEl2_lasso,1);...
               mean(WeightedSFEl2_post,1);...
               mean(WeightedSFEl2_bic_raw,1);...
               mean(WeightedSFEl2_bic_std,1);...
               mean(WeightedSFEl2_ar,1)]/VAR1_lasso,'-o')
legend({'VAR(q) Lasso';'VAR(q) Post-Lasso';...
        'VAR(q) BIC-Lasso (raw)';'VAR(q) BIC-Lasso (std)';...
        'Univ. AR(q)s'});
%%

[mean(WeightedSFEl2_ar,1);mean(WeightedSFEl2_ols,1);mean(WeightedSFEl2_lasso,1);mean(WeightedSFEl2_post,1);mean(WeightedSFEl2_bic,1)]

%%
[mean(WeightedSFEl2_ar,1);mean(WeightedSFEl2_lasso,1);mean(WeightedSFEl2_post,1);mean(WeightedSFEl2_bic,1)]
%%
100*[mean(WeightedSFEl2_lasso,1);mean(WeightedSFEl2_post,1);mean(WeightedSFEl2_bic,1)]/mean(WeightedSFEl2_lasso(:,1),1)
%%
100*mean(WeightedSFEl2_ols(:,1),1)/mean(WeightedSFEl2_lasso(:,1),1)
%% Relative Forecasting Performance (VAR(1) Lasso as benchmark)
dosave=1;
h=figure;
plot(100*[mean(WeightedSFEl2_lasso,1);mean(WeightedSFEl2_post,1);...
          mean(WeightedSFEl2_bic_raw,1)]'/mean(WeightedSFEl2_lasso(:,1)),...
          '-o','LineWidth',1);
legend({'$\mathrm{VAR}(q)$ Lasso';...
        '$\mathrm{VAR}(q)$ Post-Lasso';'$\mathrm{VAR}(q)$ BIC-Lasso'},...
        'interpreter','latex','location','east');
xlim([1,12]);
xlabel('Order of Autoregression, $q$','interpreter','latex');
ylabel('MIVWSFE in \% of $\mathrm{VAR}(1)$ Lasso','interpreter','latex');
% title('FRED-MD: Forecasting Performance','interpreter','latex');
grid on;
% => minimized at VAR(1) Post-Lasso
if dosave
    cd 'D:\Dropbox\High-Dimensional Time Series\Matlab\hdvar\application\FRED\img'
    filename=sprintf('%s_Figure_FRED-MD_Relative_MIVWSFE_lasso_post_bic', date);
    saveFigure(h,filename,'portrait',1/3)
end


