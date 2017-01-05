
function result =empirical_map(x_tr,y_tr,x_val,y_val,order,es,par_state,verbosity)
%% define logging levels
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;
%%

%empirical 

[Md_t,H]=createDesignMatrix(x_tr,es,order);

% if thre is a valida
[Md_v,H]=createDesignMatrix(x_val,es,order);
x_val=flatten_t(x_val,H);
y_val=flatten_t(y_val,H);
x_tr=flatten_t(x_tr,H);
y_tr=flatten_t(y_tr,H);
% cov function
Kov_t=Md_t*Md_t';
% estimating Hyperparameter
 msg = 'Computing Hyperparameter';

 console_output(verbosity, msg, LOG_INFO_MAJOR);
 
[tauq, sigq] =estHyperParameter(Kov_t,y_tr,par_state);


if tauq<=0.0
    
    
    tauq=1e-20;
end

l2=sigq./tauq;

%sinbular value decomposition
[xx,ss,vv]=svd(Md_t,'econ');
ss_e=diag(ss);
ss1=diag(ss_e./(ss_e.^2+l2));
w= (vv*ss1)* (xx'*y_tr);
y_pred_t=Md_t*w;

rr_t=r_2(y_pred_t,y_tr,0);
% structure of result re
result.datatype='Gauss_pro';
result.regtype='MAP';
result.order=order;
result.optimal_w=w;
result.tau=tauq;
result.sigma=sigq;
result.l2_regul=l2;
ss1=diag(ss_e./(ss_e.^2+result.l2_regul));
Hii= diag((Md_t*((vv*ss1)*xx')));
result.yPre_cv=(y_pred_t-(Hii.*y_tr))./(1-Hii);
result.r_2_cv=r_2(result.yPre_cv,y_tr,0);
result.nmrse_cv=nrmse1(result.yPre_cv,y_tr);
result.nmrse_t=nrmse1(y_pred_t,y_tr);
result.rr_t=rr_t;
% if x_val is present
%result for x_val if not skip

result.yPre_val=(Md_v*w);
result.r_2_val=r_2(result.yPre_val,y_val,0);
result.nmrse_val=nrmse1(result.yPre_val,y_val);
result.x_val=x_val;
result.y_val=y_val;
result.design_val=Md_v;
result.x_tr=x_tr;
result.y_tr=y_tr;
result.design_tr=Md_t;
result.y_pred_t=y_pred_t;
result.svd_x=xx;
result.svd_s=ss_e;
result.svd_v=vv;
try
    result.varp=result.sigma.*Hii;
catch
end






