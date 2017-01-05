


function norm = normalize_data(x)
    
    trials=size(x,2);
    data_points=size(x,1);
    x=reshape(x,size(x,1)*size(x,2),1);
    m=prctile(x,[5,95]);
    norm=(x-mean(x))/(m(2)-m(1));
    norm=reshape(norm,data_points,trials);
       
    
end

