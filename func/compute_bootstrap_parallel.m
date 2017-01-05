function[nrmse_b]=compute_bootstrap_parallel(nci,rdu,nr,processed_data)
nrmse_b=NaN(nci,1);

parfor j=1:nci
%                  u =ceil(nr*rand(nr,1));
                 ak=datasample(rdu,nr);
%                  ak=rdu(indx(j,:));
                 yP_b=bsxfun(@plus,processed_data,ak);
                 
%                 nrmse_b(j,:)=nrmse1(Result.est_param.channel(it).test(iter).delay(i).yPre_cv,yP_b(:,:)');
                
                nrmse_b(j,:)=nrmse1(processed_data,yP_b(:,:));
               
end
end