
%create kernel matrix

function [K,H] = createKernelMatrix(x1,x2,es,order)
% 

[kseq,subset]=convertES(es(1),es(2));

%if length(x1)==length(x2)
if size(x1,2)==size(x2,2)
   equalx='True';
else
    equalx='False';
end

H=kseq;
n_F=kseq+1;




if size(x1,1)==1
    n_tr=1;
    n_t=size(x1,2)-H;
    nt=n_t;
else
    n_tr=size(x1,1);
    nt=size(x1,2)-H;
    n_t=nt*(size(x1,1));

end

if strcmp( equalx,'False' )
    if size(x2,1)==1
        n_vr=1;
        n_v=size(x2,2)-H;
        nv=n_v;
    else
        n_vr=size(x2,1);
        nv=size(x2,2)-H;
        n_v=nv*(size(x2,1));
    end
end

    
    
if  strcmp( order,'inf' )  %'infinite' volterra
    X_1=zeros(n_t,n_F);
    idx=1;
    for k =1:n_tr
        if n_tr==1    % if there is only one trial 

            for i=1:nt
                X_1(idx,:)=x1(i:i+H);
                idx=idx+1;
            end
    
        else
            for i=1:nt
                X_1(idx,:)=x1(k,i:i+H);
                idx=idx+1;
            end
       end
    end
   
    if es(2)>1
        X_1=X_1(:,subset);
    end
    if strcmp( equalx,'True' )
%         K=exp(dot(X_1,X_1));
        K=exp(X_1*X_1');  

        return 
       
    end
   
    X_2=zeros(n_v,n_F);
    idx=1;
    for k =1:n_vr
    if n_vr==1    % if there is only one trial 

        for i=1:nv
            X_2(idx,:)=x2(i:i+H);
            idx=idx+1;
        end

    else
        for i=1:nv
            X_2(idx,:)=x2(k,i:i+H);
            idx=idx+1;
        end
   end
   end
   if es(2)>1
       X_2=X_2(:,subset);
   end
   if strcmp( equalx,'False' )
       K=exp(X_1*X_2');
%        K=exp(dot(X_1,X_2'));
       return 
    end
end

if strcmp( order,'inf' )==0
    
    [Md_1,H]=createDesignMatrix(x1,es,order);
    if strcmp(equalx,'True')
          K=(Md_1* Md_1');
%         K=dot(Md_1, Md_1');
    else
        
        [Md_2,H]=createDesignMatrix(x2,es,order);
%         K=dot(Md_1, Md_2');
        K=Md_1* Md_2';
   end
end



        

%if not equal do the same