
function [Result] = DDIFT_delays_analysis(varargin)

% DDIFT_delays_analysis performs the main analysis routine.
% You can call this function as follows: 
%
%  Result=DDIFT_delays_analysis(data_prepare)
%
%
%*** INPUT PARAMETERS  
%
%    data_prepare =  output of the function DDIFT_prepare.m
%
%*** OUTPUT 
%
%    Result= struct() with the following field:
%
%    .dlResult=(n.cadidates delays X 6)dlResult. The matrix is generated  
%     for each channel combination tested and for both direction 
%     of interactions (cfg.Testing=”True”).
%
%     nx1 = Training set NRMSE
%     nx2 = Cross validation NRMSE
%     nx3 = Validation setNRMSE
%     nx4 = Lower confidence interval NRMSE
%     nx5 = Upper confidence interval NRMSE
%     nx6 = Baseline NRMSE
%
%
%    .est_param=substructure with the following fields:
%
%        .Delay= Delay is a structure. It contains parameters
%                estimated for each candidate delay ? 
%
%        .param_ci= [1x2] Upper and lower confidence intervals for 
%                   the estimated delay ?
%
%        .param_final_delay= The estimated delay ?
%
%        .param_pNrmse_t=Interpolation of NRMSE training
%
%        .param_pNrmse_cv=Interpolation of NRMSE cross-validation
%
%        .param_pNrmse_lower=Interpolation of NRMSE lower confidence interval
%
%        .param_pNrmse_upper=Interpolation of NRMSE upper confidence interval
%
%
%
%    .cfg=(substructure)It keeps basic information of  the analysis used 
%         for easy plotting
%
%        .channelcombination={n.channel combination x 2}Cell array with 
%                             channel  labels that specifies  the channel 
%                             combinations  analyzed
%
%
%        .delay=Range of tested delays in ms
%
%        .iteration=It indicates how many channels pairs were analysed
%
%        .test=It indicates if one or both directions of interaction 
%              were tested.    
%%
%Quick check if the user is providing the correct input and
%using prepared  data with D^2IFT.prepare.m

if isfield(varargin{1},'DDIFT_prepare')
    
    data = varargin{1};
else
    error('\DDIFTOOL: incorrect input values, You must provide the output from DDIFT_prepare.m!');
end
%%

LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;
% LOG_DEBUG_COARSE = 3;
% LOG_DEBUG_FINE = 4;
%inizitialize a structure of results 
time_start=tic;
Result=structure_result(data,data.DDIFT_prepare.cfg);

% indices of selected channel, keep combination
[~,indices,list_ind]=channel_select(data,data.DDIFT_prepare.cfg);
Result.cfg.channelcombination=list_ind;

%for each combination take the data from the channel to analysed channels pair 
it=1;
for f_c=1:size(indices,1)
    f_ch=indices(f_c,1);
    data1=data.DDIFT_prepare.pre_ch(f_ch).channel;
    data1=data1';
    
    t_ch=indices(f_c,2);
    data2=data.DDIFT_prepare.pre_ch(t_ch).channel;
    data2=data2';
    % keep name of channels for plotting and console
    vv=list_ind{f_c,1};
    vv1=list_ind{f_c,2};
  
    
    % analysis on both  or one direction 
    iter=1;
    while 1
        if iter==1
            interaction='A-B';
            disp_inter=['direction of interaction ', vv, '-' ,vv1];
        elseif iter>1 && strcmp(data.DDIFT_prepare.cfg.testing,'True');
            interaction='B-A';
            disp_inter=['direction of interaction ', vv1, '-' ,vv];
        else
            break
        end
        
        msg =[ 'Testing channel ' ,vv,' - ' ,vv1, ' / ',disp_inter];

        console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
        maxRec = 1;
        
      
        %fit statistical model for each candidate delay
        for i=1:length(data.DDIFT_prepare.cfg.delay)
              
%            
            
            switch(interaction)
            case 'A-B'
                d_shift=data.DDIFT_prepare.cfg.delay(i);
                if d_shift>0
                    driverM=data1(1: end -d_shift,:);
                    slaveM=data2(d_shift+1: end ,:);
                else 
                    driverM=data1;
                    slaveM=data2;
                end
            case 'B-A'

                d_shift=data.DDIFT_prepare.cfg.delay(i);
                if d_shift>0
                    driverM=data2(1: end -d_shift,:);
                    slaveM=data1(d_shift+1: end,: );
                else 
                    driverM=data2;
                    slaveM=data1;
                end
            end
            
            msg =[ 'Testing delay - ' ,num2str(d_shift)];

            console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
         

            %normalize data with function normalize_data.m
            if strcmp( data.DDIFT_prepare.cfg.order,'inf')
                % noralization with zscore gave error. USE now normalize.m
%                 driverM=zscore(driverM);
%           
%              
%            
%                 slaveM=zscore(slaveM);
                 driverM=normalize_data(driverM);



                 slaveM=normalize_data(slaveM);

%                 
            else
                
                driverM=normalize_data(driverM);



                slaveM=normalize_data(slaveM);
            end

            %split data in training  and validation set

            n=max(round(size(driverM,1)/2),size(driverM,1)-data.DDIFT_prepare.cfg.numvalidate);

            x_tr=driverM(1:n,:);
            y_tr=slaveM(1:n,:);
            x_val=driverM(n+1:end,:);
            y_val=slaveM(n+1:end,:);

            % baseline_nrmse(x_t,y_t)
            %define lag
            bl=min(ntmses_b2(y_tr,x_tr,round(length(x_tr)/2)-1));
            % baseline estimation goes in dlResult in the 6 th column
            Result.dlResult.channel(it).test(iter).result(i,6)=nrmse1(y_tr,x_tr);

            x_tr=x_tr';
            y_tr=y_tr';
            x_val=x_val';
            y_val=y_val';
            
            %statistical model 
            if strcmp( data.DDIFT_prepare.cfg.order,'inf')
                 %full expansion volterra series
                 try
                    Result.est_param.channel(it).test(iter).delay(i)=empirical_GPR(y_tr,x_tr,y_val,x_val,data.DDIFT_prepare.cfg.order,data.DDIFT_prepare.cfg.es,data.DDIFT_prepare.cfg.par_state,data.DDIFT_prepare.cfg.verbosity);
                 catch
                     fprintf('\n')
                     error('DDIFTOOL error: there are not enough data points for this set of candidate delays and parameters ');
                 end

            else
                % volterra series for order 1,2,3,..
                try
                    Result.est_param.channel(it).test(iter).delay(i)=empirical_map(y_tr,x_tr,y_val,x_val,data.DDIFT_prepare.cfg.order,data.DDIFT_prepare.cfg.es,data.DDIFT_prepare.cfg.par_state,data.DDIFT_prepare.cfg.verbosity);
                catch
                    fprintf('\n')
                    error('DDIFTOOL error: there are not enough data points for this set of candidate delays and parameters ');
                end
            end
            
            % add results to the structure Result
            Result.dlResult.channel(it).test(iter).result(i,1)=Result.est_param.channel(it).test(iter).delay(i).nmrse_t;
            Result.dlResult.channel(it).test(iter).result(i,2)=Result.est_param.channel(it).test(iter).delay(i).nmrse_cv;
            Result.dlResult.channel(it).test(iter).result(i,3)=Result.est_param.channel(it).test(iter).delay(i).nmrse_val;
            
            %parametric bootstap of residuals
  
            nci=50000;
            rdu=Result.est_param.channel(it).test(iter).delay(i).y_tr-Result.est_param.channel(it).test(iter).delay(i).yPre_cv;
            nr=length(rdu);
            
%            
            processed_data=Result.est_param.channel(it).test(iter).delay(i).yPre_cv;
           
           % I checked performance with function bootstrp with a minor
           % change ( the first input is not resampled, only residual),
           % time results were the same (a bit slower).
%            opt = statset('UseParallel',true);
%            [nrmse_b]=bootstrp1(50000, @nrmse2,processed_data, rdu, 'Options', opt);
           
            if  ~isempty(data.DDIFT_prepare.cfg.par_state) 
                 %computing parfor
                  msg = 'Computing Parametric Bootstrapping of residuals';

                  console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
                
                  [nrmse_b]=compute_bootstrap_parallel(nci,rdu,nr,processed_data);
                
            else
                %for loop
                  msg = 'Computing Parametric Bootstrapping of residuals';

                 console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
              
                 [nrmse_b]=compute_bootstrap(nci,rdu,nr,processed_data);
          
           end
          
             %changed see above
%          
%            for j=1:nci
% %                 
%                  dat_samp=datasample(rdu,nr);
% %                  
%                  yP_b=bsxfun(@plus,processed_data,dat_samp);
%                  
% %                nrmse_b=compute_bootstrap(nci,it,iter,i,Result,yP_b);
%             
%            end
           
%            

                        
            pct_lower=prctile(nrmse_b,1);
            pct_upper=prctile(nrmse_b,99);
            Result.dlResult.channel(it).test(iter).result(i,4)=pct_lower;
            Result.dlResult.channel(it).test(iter).result(i,5)=pct_upper;
            
            %if NRMSE below 1 
            if Result.est_param.channel(it).test(iter).delay(i).nmrse_cv < maxRec 
                maxRec=Result.est_param.channel(it).test(iter).delay(i).nmrse_cv;
                mDelay = data.DDIFT_prepare.cfg.delay(i);
                mC=[pct_lower,pct_upper];
                mCI.ff(it).test(iter).result(1,1) = mC(1);
                mCI.ff(it).test(iter).result(1,2) = mC(2);
            else
                try 
                    if isempty(mCI.ff(it).test(iter).result(1,1)) && isempty(mCI.ff(it).test(iter).result(1,2))
                        mC=[NaN,NaN];
                    
                    
                        mCI.ff(it).test(iter).result(1,1) = mC(1);
                        mCI.ff(it).test(iter).result(1,2) = mC(2);
                    end
                catch Me 
                    mC=[NaN,NaN];
                    
                    
                    mCI.ff(it).test(iter).result(1,1) = mC(1);
                    mCI.ff(it).test(iter).result(1,2) = mC(2);
                    Result.cfg.error=Me;
             
                
                end
                
            end

        end
        
        iter=iter+1;
        if iter>2
            msg = 'Testing causual and a-causal hypothesis done';

            console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
            break
        end

     end
    
        
    
    it=it+1;
end


it=it-1;
Result.cfg.iteration=it;
%% save result later
% assume data_DDIFT already exist after DDIFT_prepare.m is called
% save result in the same subfolder of the one created by DDIFT_prepare.m name es: sub1
%

% %open directory core function and subfolder data_DDIFT
% pathDDIFT = which('DDIFT_prepare.m');
% pathDDIFT = pathDDIFT(1:end-15);
% full_n2=fullfile(pathDDIFT,'data_DDIFT\',data.DDIFT_prepare.cfg.fold_name);
% % or use cd to open the folder
% filename=[full_n2 '\' 'Result_analysis'];
% save(filename,'Result');
% save result
% save([full_n2 data.DDIFT_prepare.cfg.fold_name '\Result_analysis' ], 'Result');
%% plot REG here

for df=1:it
    
    msg = 'Preparing REG plot';

    console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
    vv=list_ind{df,1};
    vv1=list_ind{df,2};
        if strcmp(data.DDIFT_prepare.cfg.testing,'True')
            plot_it=2;
        else
            plot_it=1;
        end
        Result.cfg.test=plot_it;
        for iter=1:plot_it
%             scale=1/data.DDIFT_prepare.cfg.sampling;
            delays=(data.DDIFT_prepare.cfg.delay*data.DDIFT_prepare.cfg.sampling*data.DDIFT_prepare.cfg.subsampling*1000)';%


            t=linspace(delays(1),delays(end),delays(end)*100+1)';
         
           
            % interpolate 
            
            pNrmse_t=splineFit(t,delays,Result.dlResult.channel(df).test(iter).result(:,1),4);
            pNrmse_cv=splineFit(t,delays,Result.dlResult.channel(df).test(iter).result(:,2),4);
            pcI_lower=splineFit(t,delays,Result.dlResult.channel(df).test(iter).result(:,4),4);
            pcI_upper=splineFit(t,delays,Result.dlResult.channel(df).test(iter).result(:,5),4);

            delEst=1;
            if any(isnan(mCI.ff(df).test(iter).result))|| isempty((mCI.ff(df).test(iter).result)) %xxxx
                delEst=0;
            else    
                if (pcI_lower(delEst)> mCI.ff(df).test(iter).result(1,2))==0 
                    delEst=0;
                else
                    while (pcI_lower(delEst)> mCI.ff(df).test(iter).result(1,2))
                        delEst=delEst+1;
                        if delEst==length(pcI_lower)
                            delEst=mDelay*100;
                            break
                        end
                    end
                end
            end

            %Corner Detection
            
            if delEst>0
                p_Dvt=gradient(pNrmse_cv);
                
                p_Dvt2=gradient(p_Dvt);
                C=p_Dvt2-20*(abs(p_Dvt)).^2;
               
                [~, row_pos]=min(p_Dvt);
%                 if row_pos==1
%                     row_pos=2;
%                 end
                
              
                thr=0-0.3*std(p_Dvt);
                [~, rb]=min(p_Dvt(row_pos:end)<thr);
                %if row_pos last element error it will crash:
                if length(p_Dvt)==row_pos
                    rb1=rb+row_pos-1;
                else
                    rb1=rb+row_pos;
                end
                [~, C_idx]=max(C(row_pos:rb1));%
                if length(p_Dvt)==row_pos
                    C_idx1=C_idx+row_pos-1;
                else
                    C_idx1=C_idx+row_pos;
                end
                
                
                
                thr2=C(C_idx1)-std(C)*0.2;
                tmp=C(C_idx1:rb1)<thr2;
                if all(tmp)==0
                 
                    [~, aka]=max(tmp);
                    C_idx2=(aka)+C_idx1- 1;%
                else
                    C_idx2=rb1;
                end
                    
               
%                 C_idx2=rb1;
                
              
                C_pos=t(C_idx2);
                HRC2=C_pos;
                HRC=t(row_pos);
                meanHRC=(mean([HRC,HRC2]));  


            else
                HRC=0;
                HRC2=0;
                meanHRC=0;
            end
        
            Result.est_param.channel(df).test(iter).param_ci=[HRC,HRC2];
            Result.est_param.channel(df).test(iter).param_final_delay=meanHRC;
            Result.est_param.channel(df).test(iter).param_pNrmse_t=pNrmse_t;
            Result.est_param.channel(df).test(iter).param_pNrmse_cv=pNrmse_cv;
            Result.est_param.channel(df).test(iter).param_pcI_lower=pcI_lower;
            Result.est_param.channel(df).test(iter).param_pcI_upper=pcI_upper;
            %% save result
            
            % assume data_DDIFT already exist after DDIFT_prepare.m is called
            % save result in the same subfolder of the one created by DDIFT_prepare.m name es: sub1
            %
            %open directory core function and subfolder data_DDIFT
            pathDDIFT = which('DDIFT_prepare.m');
            pathDDIFT = pathDDIFT(1:end-15);
            full_n2=fullfile(pathDDIFT,'data_DDIFT\',data.DDIFT_prepare.cfg.fold_name);
            % or use cd to open the folder
            filename=[full_n2 '\' 'Result_analysis'];
            save(filename,'Result');
                       
            
            % XXX save parameter Hrc hici loci and mean
            % interpolation and data for future plotting is done
            % is user indicated display, do it here.
            if  strcmp(data.DDIFT_prepare.cfg.display,'True') 
                
                if strcmp(data.DDIFT_prepare.cfg.testing,'True')==1 && iter==1

                    tit=['Results channels ',vv,'->' ,vv1 ] ;
                    pl=['REG',vv,vv1];
                elseif strcmp(data.DDIFT_prepare.cfg.testing,'True')==1 && iter==2
                    tit=['Results channels ', vv1,'->' ,vv ] ;
                    pl=['REG',vv1,vv];
                else
                    tit=['Results channels ',vv,'->' ,vv1] ;
                    pl=['REG',vv,vv1];
                end
    
                hFig = figure(df);    
                subplot(1,2,iter);
                title(tit,'Color', 'r')
                set(hFig, 'Position', [10 10 10 10 ])
                hold on 
                plot(t,pNrmse_t,'--','linewidth',2,'Color',[0 0 0]+0.05);
                plot(t,pNrmse_cv,'linewidth',2,'Color',[0 0 0]+0.05);
                plot(t,pcI_lower,'linewidth',1,'Color',[0 0 0]+0.05);
                plot(t,pcI_upper,'linewidth',1,'Color',[0 0 0]+0.05);

                plot(delays,Result.dlResult.channel(df).test(iter).result(:,2),'o','linewidth',1,'Color',[0 0 0]+0.20);
                fill([t ;flipud(t)], [pcI_lower; flipud(pcI_upper)],'k');
                alpha(0.25);
                ylim([0 1]);
                xlim([t(1),t(end)]);
 
                vline(HRC,'k');
                vline(HRC2,'k');
                vline(meanHRC,'k');
                % do not crash for old matlab, where round take only one
                % input
                try
                    text(HRC2+0.01,0.12,['ci_u', ' -- ' , num2str(round(HRC2,1)),' ms']);
                    text(HRC2+0.01,0.15,['Delay', ' -- ' , num2str(round(meanHRC,1)),' ms'],'Color','red','FontSize',10);
                    text(HRC2+0.01,0.18,['ci_l', ' -- ' , num2str(round(HRC,1)),' ms']);
                catch Me
                    text(HRC2+0.01,0.12,['ci_u', ' -- ' , num2str(HRC2),' ms']);
                    text(HRC2+0.01,0.15,['Delay', ' -- ' , num2str(meanHRC),' ms'],'Color','red','FontSize',10);
                    text(HRC2+0.01,0.18,['ci_l', ' -- ' , num2str(HRC),' ms']);
                end
                grid('on');
                xlabel('Delay');
                ylabel('NRMSE');
                hold off
                % save figure in the folder
                name=pl;
                saveas(hFig,[full_n2  '\' name],'fig')
            end
        end
%%
 time_end=toc(time_start);
 msg = sprintf( ...
    'DDIFTOOL analysis ended: %s \n analysis took %.0f MINUTES (%.0f SECONDS)', ...
    datestr(now), time_end/60, time_end);
 console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
%stop parallel computing	
delete(data.DDIFT_prepare.cfg.par_state);

end








    