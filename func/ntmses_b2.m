
function base_n = ntmses_b2(x,y,lag)


%Computes the normalized root mean squared error of 
%time series x and y, possibly with delay shift 'lag'.



if lag==0
    base_n=nrmse1(x,y)
else
    
   le=2*lag+1;
   r=zeros(le,1);
   for i =-lag:lag
       
      
       if i>0
           r(i+lag+1,:)=nrmse1(x(1:end-abs(i)),y(abs(i)+1:end));
       elseif i==0
           r(i+lag+1,:)=nrmse1(x,y);
       
       else
           r(i+lag+1,:)=nrmse1(y(1:end- abs(i)),x(abs(i)+1:end));
       end
   end
end
   base_n=r;
