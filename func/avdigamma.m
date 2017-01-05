function[avg] = avdigamma(sample,vect)

N=length(sample);
atria = nn_prepare( sample,'maximum');
avg=0.0;
ncount = zeros(length(sample),1);
for i=1:N
    dist=vect(i);
    
    [count, neighbors] = range_search(sample,atria,i,dist-1e-15,0);
    
    ncount =count;
   
    avg=avg+psi(count)/N;
    
end
end