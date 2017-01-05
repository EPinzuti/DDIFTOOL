%bspline

 function bs = bspline(x,t,i,r)


    bs=zeros(length(x),1);

    for idx=1:length(x)
        bs(idx,:)= CDB(x(idx),t,i,r);
    end
