
function yt = splineFit(x,t,y,o)



S=spline_t(x,t,o);

S=S';
idx=[];
k=1;
for i=1:length(x)
    if abs(x(i) -t(k))<1e-10
        idx=[idx;i];
    else
        continue
    end
    k=k+1;
  
end

y2= sf(x,t,y,idx);
[x,s,v]=svd(S,'econ');

s_e=diag(s);
s1=diag(s_e./(s_e.^2));
w=(v*s1)*(x'*y2);
yt=S*w;