

function y = flatten_t(y,H)

if (size(y,1))==1
    y=y(H+1:end);
    y=y';
else
    tm=(y(:,H+1:end));
    tm=tm';
    y=reshape(tm,size(tm,1)*size(tm,2),1);
end
    
 