

function [Md_t, H] = create_big_DesignMatrix(x_tr,knots,order,subset)
%DESIGN MATRIX [for n-th order Homogeneous Polynomial System]

% if 

% else

H=knots;
nF=knots+1;
n_t=size(x_tr,2)-H;
H1_F=zeros(n_t,nF);
for i=1:n_t
    
    H1_F(i,:)=x_tr(i:i+H);
end

H1_F=H1_F(:,subset);
nF=size(H1_F,2);
if order==1
    Md_t=[ones(n_t,1),H1_F];
    return 
end

idcs=ones(nF,order);
for k=1:order
    for n=1:nF
        if k-1==0
            idcs(n,k)=((nchoosek(1,k)));
        else
            
            idcs(n,k)=nchoosek(n-1+k-1,k-1);
        end
    end
end

n_M=zeros(1,order);
for ty=1:order
    n_M(:,ty)=sum(idcs(:,ty));
end
n_M=sum(n_M)+1;
Md_t=ones(n_t,n_M);
ftable={};
for i=1:2
   
    ftable.f{i}=zeros(sum(idcs(:,i)),i);
    idx=1;
    for j=1:nF
       
        ftable.f{i}(idx:idx-1+idcs(j,i),1)=j;
        idx=idx+idcs(j,i);
    end
end
idx=1;
for j=1:nF
   

    ftable.f{2}(idx:idx-1+idcs(j,2),2)=ftable.f{1}(1:idcs(j,2),:);
    idx=idx+idcs(j,2);
  
end
   
idx=2;
for i=1:2
    
    didx=sum(idcs(:,i));
    m_remap=reshape(H1_F(:,ftable.f{i}),[length(H1_F(:,1)),length(1:didx),length(1:i)]);
    %Md_t(:,idx:idx-1+didx)=prod(H1_F(:,ftable.f{i}),3);
    Md_t(:,idx:idx-1+didx)=prod(m_remap,3);
    idx=idx+didx;
end
%for order grater then 2

for i=3:order
    

    didx=sum(idcs(:,i));
    pidx=sum(idcs(:,i-1));
    ridx=[]; 
    cidx=[];
    for j=1:nF
       se=linspace(1,idcs(j,i),idcs(j,i));
       ridx=[ridx,se];
       tee=repmat(j,1,idcs(j,i));
       cidx=[cidx,tee];
       
    end
    stp=Md_t(:,idx-pidx:idx-1);
    Md_t(:,idx:idx-1+didx)=H1_F(:,cidx).* stp(:,ridx);
    idx=idx+didx;
   
   
    
end

