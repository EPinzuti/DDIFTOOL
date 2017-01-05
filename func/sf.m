

%sf function
function sp = sf(x2,t,y,idx)
S=spline_t(x2,t,2);

S=S';

[x,s,v]=svd(S(idx,:),'econ');

s_e=diag(s);
s1=diag(s_e./(s_e.^2));
try
    w=(v*s1)*(x'*y);
catch Me
    fprintf('\n')
    error(' error: there is a bug that have to be fixed!: Please change slightly delay candidate  or subsampling');
end
% w=dot(dot(v',s1),dot(x',y))
sp=S*w;

