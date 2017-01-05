function[paral_state poolobj]= set_parallel(cfg)
%the function check if is possible to start parallel computing
%Part of the code is taken from set_parallel (TRENTOOL)
ver = version('-release');
if str2double(ver(1:4)) >= 2014
    ver = 'newer';
    
elseif str2double(ver(1:4)) < 2013
    ver = 'older';   
else
end    
if  ft_hastoolbox('DCT');
            
    paral_state=1;
        
        switch ver
            case {'newer'}
                parallelConfig = parcluster(parallel.defaultClusterProfile);
            case {'2013a' '2013b' 'older'}
                parallelConfig = findResource('scheduler','configuration',defaultParallelConfig);
        end
        workers = parallelConfig.NumWorkers;  
        ver = version('-release');
        if str2double(ver(1:4)) >= 2014 
            poolobj = gcp('nocreate'); 
            if isempty(poolobj)
                poolsize = 0;
                parpool(workers,'IdleTimeout', Inf)% If no pool, do  create new one.
                poolobj = gcp('nocreate');
            else
                poolsize = poolobj.NumWorkers;
            end
        else
             poolobj = gcp('nocreate'); 
            if isempty(poolobj)
                poolsize = 0;
%                 matlabpool(workers,'IdleTimeout', Inf)% If no pool, do create new one.
            else
                poolsize = poolobj.NumWorkers;
            end
        end
else
    paral_state=0;
end
end
        
        
        
        

