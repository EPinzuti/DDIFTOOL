
% create design matrix fro volterra series
function [B_design,H1] = createDesignMatrix(x,es,order)


%check if a es is a list major of twoa

%if not

[kseq,subset]=convertES(es(1),es(2));
%check if there is one or more trials, for now we deal with one trial

%if length(size(x))==2
if size(x,1)==1
    [B_design, H1]=create_big_DesignMatrix(x,kseq,order,subset);
else 
    n_tr=size(x,1);
    B_design=[];
    %B_design=cell(1,n_tr);
    for k=1:n_tr
        [temp, H1]=create_big_DesignMatrix(x(k,:),kseq,order,subset);
        B_design=[B_design;temp];
        %B_design{k}=temp;
    end
    %B_design = cell2mat(B_design);
        
        
    
end


end


