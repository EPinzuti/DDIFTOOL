


function mi =mutual_inf(x,y,k)
% mutual information based on NON-PArametric Entropy Estimation (Ver Steeg, 2000)
%nn_search and nn_prepare are OpenTStool functions to estimate k-nearest neighboors 
%efficiently

% k=3
%  x=y(1:end-1,:);
%  y=y(1+1:end,:);
int=1e-10;
x1=zeros(length(x),1);
y1=zeros(length(y),1);
for p=1:length(x)
    
    x1(p,:)=x(p)+int.*rand(1);
    y1(p,:)=y(p)+int.*rand(1);
end

points=[x1,y1];
atria = nn_prepare( points,'maximum');

[index, distance] = nn_search( points, atria, 1:length(x), k,0);

dvec=distance(:,3);
a=avdigamma(x1,dvec);
b=avdigamma(y1,dvec);
c=psi(k);
d=psi(length(x1));
mi=(-a-b+c+d)/log(2);




    
