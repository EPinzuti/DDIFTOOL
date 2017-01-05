
function det_subsampling(x)
 %Determines by mutual information a feasible subsampling of a given 
 % time series x.
 
 wind=50;
 x=driverM
 n=min(n,length(x));
 mia=zeros(wind,1);
 for i= 1:wind
     mia(i,:)=mi(abs(x(1:end-i))  ,abs(x(i+1:end)));
     
 end
     
 
 %use function amutual from tstool
 
 
 