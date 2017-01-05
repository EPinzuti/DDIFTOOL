
% Computes BCa  NRMSE confidence intervals for gaussian process object gpo,     on a subset of n data points, using B bootstrap sample sets created by
% non-parametric resampling of LOO-Predictive distribution.

function yp =evaluate_at(result, tv_t,idx_t,coeff_t)
%gpo is a structure with result from delay analysis
%idx indices of data point which are used for evaluation
% determined regression coefficients
if nargin < 1
    error('neco: Please provide ''result'' as input!');
end
gpo=result;
if  nargin == 1
    tv_t='t';
    idx_t='False';
    coeff_t='False';
elseif nargin==2
    idx_t='False';
    coeff_t='False';
elseif nargin==3
    coeff_t='False';

else
    
end

model_type=gpo.delay.regtype;
if strcmp(model_type,'MAP')
    if strcmp(idx_t,'False')
        if strcmp(tv_t,'v')
            Md_b=gpo.delay.design_val;
        elseif strcmp(tv_t,'t')
            
            Md_b=gpo.delay.design_tr;
        else
        end
    else
        idx=idx_t(1);
        if strcmp(tv_t,'v')
            Md_b=gpo.delay.design_val(idx,:);
        elseif strcmp(tv_t,'t')
                       
            Md_b=gpo.delay.design_tr(idx,:);
        else
        end
    end
    if strcmp(coeff_t,'False')
        wb=gpo.delay.optimal_w;
    else
        wb=coeff_t;
    end
    yp=Md_b*wb;
    
else
end

if strcmp(model_type,'GPR')
    
    if strcmp(idx_t,'False')
        if strcmp(tv_t,'v')
            K_b=gpo.delay.kerneldesign_val;
        elseif strcmp(tv_t,'t')
            
            K_b=gpo.delay.kerneldesign_tr;
        else
        end
    else
        if strcmp(tv_t,'v')
            K_b=sker(gpo.delay.kerneldesign_val,idx);
        elseif strcmp(tv_t,'t')
            
            K_b=sker(gpo.delay.kerneldesign_tr,idx);
        else
        end
    end
    if strcmp(coeff_t,'False')
        if strcmp(idx_t,'False')
            ab=gpo.delay.optimal_hyp;
        else
            ab=sker(gpo.delay.optimal_hyp,idx);
        end
    else
        ab=coeff;
    end
    
    yp=K_b*ab;
    
end
    
% XXX check idx is correct need to put idx start- end    
    


        

    
    
    