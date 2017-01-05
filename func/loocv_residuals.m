
%loo-cv residuals
function lloo= loocv_residuals(K,y,lam)

I=eye(size(K,1));



try
    [R p] = chol(K+lam.*I);
    Pinv= R\(R'\I);
catch Me
    [K,te]=simpleNearpd(K);
    
%     K=nearestSPD(K);
    R  = chol(K+lam.*I);
    %R = chol(K.*tau0+I);
    Pinv= R\(R'\I);
    
end
%Pinv= R\(R'\I);
%if error add positive value to eigenvalue

pii=diag(Pinv);
mu=(Pinv*y)./pii;
lloo=sum(mu.^2)./length(y);
