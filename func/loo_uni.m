
%function loo_uni_
function [logpi1 ]= loo_uni(tau0,K,y)

%tau0=18.25992299
n=size(K,1);
I=eye(n);

% [R p] = chol(K.*tau0+I);   
 %[R p] = chol(bsxfun(@plus,(K.*tau0),I));
% if p==0
%    
%     Pinv= R\(R'\I);
% else 
%     
%     while p ~=0
%         [K,te]=simpleNearpd(K);
%         %[R p]  = chol(bsxfun(@plus,(K.*tau0),I));
%         [R p]= chol(K.*tau0+I);
%         Pinv= R\(R'\I);
%     end
% end
%Pinv= R\(R'\I);

try
    
    [R,p] = chol(K.*tau0+I);  
    Pinv= R\(R'\I);
     %Pinv=(chol(K.*tau0+I))\((chol(K.*tau0)+I)'\I);

catch Me
    
    [K,te]=simpleNearpd(K);
    %K=nearestSPD(K);
   
    [R,p]= chol(K.*tau0+I);
    Pinv= R\(R'\I);
   
    %Pinv=(chol(K.*tau0+I))\((chol(K.*tau0+I))'\I);

end


%if error add positive value to eigenvalue

pii=diag(Pinv);
q_n2=(Pinv*y).^2;
sig=sum(q_n2./pii)/n; 
logpi=(1./(2*n*sig))*(n*sig)-(1./(2*n))*sum(log(pii));
logpi1=(logpi+0.5*log(2*pi*sig));


 

