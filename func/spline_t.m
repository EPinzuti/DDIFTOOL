

%spline for regression graph

function spline1=spline_t(x,t,r)
    
  
    
   t_ext=t;
   
   init=repmat(t_ext(1),(r-1),1);
   final=repmat(t_ext(end),(r-1),1);
   t_ext =[init; t_ext;final];
   
   n = length(t_ext)-r;
   spline1 = zeros(n,length(x));
   for i=1:n
       spline1(i,:)= bspline(x,t_ext,i,r);
  
   
   end
  
 end