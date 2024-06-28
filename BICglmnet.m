function bic = BICglmnet(X,Y,glmnet_options)
fit=glmnet(X,Y,'gaussian',glmnet_options);
RSS=sum((Y-fit.a0'-X*fit.beta).^2,1)';
n=size(X,1);
BIC=n*log(RSS)+log(n)*fit.df;
[~,id_min]=min(BIC);
lambda=fit.lambda(id_min);
beta=fit.beta(:,id_min);
if glmnet_options.intr
    a0=fit.a0(id_min);
else
    a0=[];
end
bic.a0=a0;
bic.beta=beta;
bic.lambda=lambda;

end