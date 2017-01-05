
function r = r_2(x,y,lag)

%Computes the squared empirical correlation coefficient of 
   % time series x and y, possibly with delay shift 'lag'.
   
    trials=size(x,2);
    data_points=size(x,1);
    x=reshape(x,size(x,1)*size(x,2),1);
    function  nm = xq(x,y)
    nm=(sum((x-mean(x)).*(y-mean(y)))/length(x)).^2./(var(x)*var(y));
    end

if lag==0
     r= xq(x,y);
else
    le=2*lag+1;
    r=zeros(le,1);
    
    for i =-lag:lag
       
      
       if i>0
           r(i+lag+1,:)=xq(x(1:end-abs(i)),y(abs(i)+1:end));
       elseif i==0
           r(i+lag+1,:)=xq(x,y);
       
       else
           r(i+lag+1,:)=xq(y(1:end- abs(i)),x(abs(i)+1:end));
       end
   end
end
end
