
function norm_nrmse = nrmse1(x,y)
%Computes the normalized root mean squared error of time series x and y

y1=reshape(y,size(y,2)*size(y,1),1);
variance=var(y1);

norm_nrmse=sqrt(sum((sum((x-y).^2)))./(length(x))./variance);%;

end