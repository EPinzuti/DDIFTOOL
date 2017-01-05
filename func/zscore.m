

function s =zscore(x)

x1=reshape(x,size(x,2)*size(x,1),1);
s=(x1-mean(x1))./std(x1);