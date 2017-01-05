

function mu_i_cv =computecv_estimator(result,idx_t,X)

model_type=gpo.delay.regtype;
if strcmp(model_type,'MAP')
    if strcmp(idx_t,'False')
                  
        Md_b=gpo.delay.design_tr;
        yb=gpo.delay.y_tr;
        
    else
        Md_b=gpo.delay.design_tr(idx_t,:);
        yb=gpo.delay.y_tr(idx_t);
            
    end
    if strcmp(X,'False')
        [xx,ss,vv]=svd(Md_b,'econ');
        ss_e=diag(ss);
        ss1=diag(ss_e./(ss_e.^2+gpo.delay.l2_regul));
    else
        xx=X;
        ss1=X;
        vv=X;
    end
    
    wb= (vv*ss1)* (xx'*yb);
    yP=Md_b*wb;
    Hii= diag((Md_b*((vv*ss1)*xx')));
    mu_i_cv=(yP-(Hii.*yb))./(1-Hii);

elseif strcmp(model_type,'GPR')
    
     if strcmp(idx_t,'False')
         
         K_b=gpo.delay.K_tr;
         yb=gpo.delay.y_tr;
         
     else
         K_b=sker(gpo.delay.K_tr,idx_t);
         yb=sker(gpo.delay.y_tr,idx_t);
     end
     I=eye(length(yb));    
     if strcmp(X,'False')
         
         [R ,p] = chol(K_b.*gpo.delay.tau+gpo.delay.sigma.*I);    
         if p==0
         

            Pinv= R\(R'\I);
         else 
            [K_b,teb]=simpleNearpd(K_b);

            %R = chol(bsxfun(@plus,K,l2).*I);
            R = chol(K_b.*gpo.delay.tau+gpo.delay.sigma.*I);
            Pinv= R\(R'\I);

         end
     else
         Pinv=X;
     end
     
     Kii=diag(Pinv);
     mu_i_cv=yb-(Pinv*yb)./Kii;
     
end

         
     

      
