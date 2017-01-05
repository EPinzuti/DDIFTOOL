
%simpleNearpd
function [K_t,t]= simpleNearpd(K_t)
% this function perform a simple correction to non-positive-definity of a square matrix
%adding small value to the diagonal until eigenvalues are positive

%quick check if there are inf or nan values, it is not use later
B = (K_t + K_t')/2;
try
    [~,~,~] = svd(B);
catch ME
    if (strcmp(ME.identifier,'MATLAB:svd:matrixWithNaNInf'))
        fprintf('\n')
        error(['DDIFTOOL error: covariance matrix can not contain inf or nan values. It is not possible to correct it to positive-definite. Make sure your'...
              'data input are correct and your are using enough data-points for analysis.']);
    end
        
end
     



n=size(K_t,1);
I=eye(n);
t=0;
%count=0;
%err_count=0;
%p=3;
[rss, p]=chol(K_t);
%if already positive-definive do not add values
if p==0
    return 
end
while p ~=0%count==err_count
    
   
     [rss, p]=chol(K_t);
    
        % is possible to use p as integer to see if it goes to an error or
        % not (0 no-error )
    %[Rss p]=chol(K_t);
     %if p==0%~isempty(Rss)
       %  break
         
        
     %else
         
    
   
    t= t+ 1e-12;
     %end
    
    %catch Me
        %err_count=err_count+1
   
    K_t=K_t+I.*1e-12;
%     %count=count+1;
%
% while 1
%     try
%         chol(K_t);
%         return
%     
%     catch Me
%         t= t+ 1e-12;
%     end
%     K_t=K_t+I.*1e-12;
 end


 
        
        
   
    


end
        
