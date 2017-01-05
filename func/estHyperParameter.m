
function [tau,sig,resl2]= estHyperParameter(K,y,par_state)
%estimate hyperparameter 
%Leave-One-Out Cross-Validation Approximation of Predictive Distribution  
%for finding hyper parameters. (paper " Predictive approaches for choosing hyperparameters in gaussian processes, Sundarajan et al. 2001) 
% K: is the kernel matrix
% y: predicted data 


rang_=min(600,length(y));
idx=linspace(1,rang_,rang_);

%K=sker(K,idx);

K=simpleNearpd(sker(K,idx));
% K=nearestSPD(sker(K,idx));

y=y(1:size(idx,2));

func =@(w)loo_uni(w,K,y);

%options = optimoptions(@fminunc,'Algorithm','quasi-newton');

%options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','StepTolerance',1e-6,'MaxFunctionEvaluations',15000,'ConstraintTolerance',1e-6);

%      minConf_TMP(func,1,1e-30,Inf);          %fmincon(func,1,[],[],[],[],1e-15,Inf,[],options)
%resl2= fminbnd(func,1e-15,20)
 %fminuc(func,1)
%resl2=fminsearchbnd(func,1,1e-15,inf)

% OPTIONS = optimset('TolX',1e-15,'TolFun',1e-15,'MaxFunEvals',10000,'MaxIter',1000);
% options = optimoptions('fmincon');
% options.Algorithm='interior-point';
% options.MaxFunctionEvaluations=15000;
% options.StepTolerance=1e-15;
% options.ConstraintTolerance=1e-15;
% 
options = optimoptions('fmincon','Display','off');
if  ~isempty(par_state) 
    options = optimoptions(options,'UseParallel',true);
else
end

% options=optimset;
% options = optimset(options,'Display', 'off') ;
% options = optimset(options,'UseParallel','always') ;

resl2 = fmincon(func,1,[],[],[],[],1e-15,inf,[],options);


%opts    = struct( 'x0',1);
% option=optimoptions(@fminunc,'Algorithm','quasi-newton','TolX',1.0000e-15);
% resl2=fminunc(func,1,option);
% opts    = struct( 'x0',1, 'pgtol', 1e-15 );
% opts.maxIts=15000;
 %resl2= lbfgsb( func, l, u );%[x,f,info] 
%did not work
%X = bfgs(1, [1 1], 1e-7, 1e-7, [1e-15 1e-15], 100, 'func')
% tau0=resl2;

tau0=resl2;
n=size(K,1);
I=eye(n);
R=chol(K.*tau0+I);
% R=chol(bsxfun(@plus,(K.*tau0),I));


% Pinv= R\(R'\I);

Pinv=solve_chol(R,I);

% Pinv=(chol(K.*tau0+I))\((chol(K.*tau0+I))'\I);

pii=diag(Pinv);
q_n2=(Pinv*y).^2;
sig=sum(q_n2./pii)./n;
tau=tau0.*sig;


