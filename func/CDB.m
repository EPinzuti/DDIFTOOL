function B_r = CDB(x,t,i,r)


    B_r=0;
    if t(i) == t(i+r)
        B_r;
        return 
    end

    if (t(i) > x) 
        B_r;
        return
    elseif  (t(i+r) < x)
        B_r;
        return 
        
    end
    if r==1
        if (x >= t(i)) && (x < t(i+1)) 
            B_r=1;
            return 
        elseif  x == t(length(t)-1)
            
            B_r=1;
            return 
            
        end
    end
    
   if (t(i+r-1)> t(i))
       B_r = B_r + (x-t(i))/(t(i+r-1) - t(i))* CDB(x,t,i,r-1) ;
   end
   if  (t(i+r) > t(i+1))
       B_r = B_r + (t(i+r)-x)./(t(i+r)-t(i+1))*CDB(x,t,i+1,r-1);
   end