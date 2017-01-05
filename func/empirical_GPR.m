function result = empirical_GPR(x_tr,y_tr,x_val,y_val,o,es,par_state,verbosity)
% empirical gpr infinite kernel regression
%% define logging levels
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;
%%

[K_tr, H]=createKernelMatrix(x_tr,x_tr,es,o);

%if x_v is not false
[K_vr]=createKernelMatrix(x_val,x_tr,es,o);
[K_vv, H]=createKernelMatrix(x_val,x_val,es,o);

x_val=flatten_t(x_val,H);
y_val=flatten_t(y_val,H);
x_tr=flatten_t(x_tr,H);
y_tr=flatten_t(y_tr,H);

n=length(x_tr);
I=eye(n);
msg = 'Computing Hyperparameter';

console_output(verbosity, msg, LOG_INFO_MAJOR);
[tauq, sigq,jj] =estHyperParameter(K_tr,y_tr,par_state);
if tauq<=0
    tauq=1e-20;
end

l2=sigq./tauq;

%[R ,p] = chol(bsxfun(@plus,K_tr,l2).*I);      
[R ,p] = chol(K_tr+l2.*I);    
if p==0
%     tic 
%     Pinv= R\(R'\I);
%     toc
%     tic 
    Pinv=solve_chol(R,I);
%     toc
else 
    [K_tr,te]=simpleNearpd(K_tr);
   
    %R = chol(bsxfun(@plus,K,l2).*I);
    R = chol(K_tr+l2.*I);
%     tic
%     Pinv= R\(R'\I);
%     toc
%     tic 
    Pinv=solve_chol(R,I);
%     toc
end

matx=Pinv*y_tr;
Kii=diag(Pinv);
mu_cv=y_tr- (matx./Kii);  % predictive mean pag 135 equation 5.12


y_pred_t=K_tr*matx;
rr_t=r_2(y_pred_t,y_tr,0);


% n-folds cross validation for now use this instead of loo-cv  since is
% less prone to overfitting, might be good to do both and see

result.datatype='Gauss_pro';
result.regtype='GPR';
result.order=o;
result.optimal_hyp=matx;
result.tau=tauq;
result.sigma=sigq;
result.l2_regul=l2;
result.r2_t=rr_t;

result.yPre_cv=mu_cv;
result.r_2_cv=r_2(result.yPre_cv,y_tr,0);
result.nmrse_cv=nrmse1(mu_cv,y_tr);
result.nmrse_t=nrmse1(y_pred_t,y_tr);

% if x_val is present
%result for x_val if not skip

result.yPre_val=(K_vr*matx);
result.r_2_val=r_2(result.yPre_val,y_val,0);
result.nmrse_val=nrmse1(result.yPre_val,y_val);
result.x_val=x_val;
result.y_val=y_val;


result.kerneldesign_val=K_vr;
result.kerneldesign_val_val=K_vv;
result.x_tr=x_tr;
result.y_tr=y_tr;
result.jj=jj;

result.y_pred_t=y_pred_t;
result.kerneldesign_tr=K_tr;

result.varp_cv=1./Kii.*result.tau; % predictive variance 117
result.varp=result.tau.*diag(K_tr-(K_tr*Pinv)*K_tr);

end




