function [ss, mi]= mutual_subsampling(y,n,k)
%call mutual_inf function to estimate mutual information 
mi=zeros(n,1);

for j= 1:n
   
    mi(j,:)=mutual_inf(y(1:end-j,:), y(j+1:end,:), k);
end
% compute optimal subsampling ( mutual information less than 1 ) ss is not used
%later. Only first minimum is taken into account
delta=zeros(n,1);
[a kl]=min(mi);
[b jl]=max(mi);
for i=1:n
    
    delta(i,:)=mi(i)-a;
end
mm=max(max(mi));
min_t=min(min(mi));
thr=0.1*(mm-min_t);
[kk ss]=max(delta<thr);
if mi(ss)>1
    if max(mi<1)
        [kk ss]=max(mi<1);
    else
        ss=ss+1;
    end
else
    ss=ss+1;
end
