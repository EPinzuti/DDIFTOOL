
%embedding space estimation sse loo-cv

function [d0,t0, gcvresid] = embed_cv2(x,y,esrange,par_state,verbosity)
%  y=data1(:,:);
%  x=data2(:,:);
%  esrange=esRange;
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;

%esrange need to be a cell array

xc=normalize_data(x');
xc=xc';
yc=normalize_data(y');
yc=yc';
if size(x,2)==1
    n=length(x);
else
    n=size(x,2);
end

dr=esrange{1};

tr=esrange{2};

gcvresid=zeros(length(dr),length(tr));
idxd=1;
for d1=1: length(dr)
    d=dr(d1);
    idxt=1;
    for t1=1:length(tr)
        t=tr(t1);
        msg =[ 'Testing dimension - ' ,num2str(d), ' and ', 'Testing tau - ' ,num2str(t) ];

        console_output(verbosity, msg, LOG_INFO_MAJOR);
        if d*t>n
            gcvresid(idxd,idxt)=inf;
            
        else 
            
            [K_t,H]=createKernelMatrix(xc,xc,[d t],'inf');
            y_tr=flatten_t(yc,H);
            func=@(lam1)loocv_residuals(K_t,y_tr,lam1);
            if n<= 500
                options = optimoptions('fmincon','Display','off');
                if  ~isempty(par_state) 
                    options = optimoptions(options,'UseParallel',true);
                else
                end
                %option=optimoptions(@fminunc,'Algorithm','quasi-newton','TolX',1.0000e-12);
                resl2=fmincon(func,1,[],[],[],[],1e-12,inf,[],options);
                lam=resl2;
            else
                lam=1.0000e-12;
            end
            gcvresid(idxd,idxt)=func(lam);
            
           
        end
        idxt=idxt+1;
    end
    idxd=idxd+1;
end

[dd, mm]=min(gcvresid(:));

doi=floor(mm./length(tr));
if doi<1
    doi=1;
end
d0=dr(doi);
t01=mm-doi.*length(tr);
if t01<1
    t01=1;
end

t0=tr(t01);  %not correct
        