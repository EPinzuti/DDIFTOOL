function hFig=REG_plot(varargin)
%%
% REG_plot allows to plot REG at later stage of analysis, once
% DDIFT_delays_analysis.m is done.(It will be useful to plot specific channel 
% pairs and subjects when the toolbox will deal with group analysis)
%
% I can be used as follow:
%
% REG_plot(data_prepare, Result)
%
% If you want to plot only a specific channel pairs change
% Result.cfg.channelcombination and provide label of channel pairs analysed
% that you want to plot
% 
%
% INPUT
% 
%  data= output of DDIFT_prepare.m
%  Reesult= output of DDIFT_delays_analysis.m
%
%
% OUTPUT
% 
% REG plot of anlysed channels
%%
%Quick check if the user is providing the correct input and

if length(varargin)==1
    error('\DDIFTOOL: incorrect number of inputs. You must provide data (from DDIFT_prepare.m) and Result !');
    
else
end
if isfield(varargin{1},'DDIFT_prepare') && isfield(varargin{2},'dlResult')
    
    data = varargin{1};
    Result=varargin{2};
else
    error('\DDIFTOOL: incorrect input.You must provide the output from DIFT_prepare.m and Result from DIFT_delays_analysis.m !');
end
%%

LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;


for df=1:Result.cfg.iteration
    msg = 'Preparing REG plot';

    console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);
    vv=Result.cfg.channelcombination{df,1};
    vv1=Result.cfg.channelcombination{df,2};
    if strcmp(data.DDIFT_prepare.cfg.testing,'True')
        plot_it=2;
    else
        plot_it=1;
    end
    for iter=1:plot_it
        delays=(data.DDIFT_prepare.cfg.delay*data.DDIFT_prepare.cfg.sampling*data.DDIFT_prepare.cfg.subsampling*1000)';

        t=linspace(delays(1),delays(end),delays(end)*100+1)';

        
        pNrmse_t=Result.est_param.channel(df).test(iter).param_pNrmse_t;
        pNrmse_cv=Result.est_param.channel(df).test(iter).param_pNrmse_cv;
        pcI_lower=  Result.est_param.channel(df).test(iter).param_pcI_lower;
        pcI_upper=Result.est_param.channel(df).test(iter).param_pcI_upper;
        HR=Result.est_param.channel(df).test(iter).param_ci;
        HRC=HR(1);
        HRC2=HR(2);
        meanHRC=Result.est_param.channel(df).test(iter).param_final_delay;
        
        if  strcmp(data.DDIFT_prepare.cfg.display,'True') 
                
                if strcmp(data.DDIFT_prepare.cfg.testing,'True')==1 && iter==1

                    tit=['Results channels ',vv,'->' ,vv1 ] ;
                elseif strcmp(data.DDIFT_prepare.cfg.testing,'True')==1 && iter==2
                    tit=['Results channels ', vv1,'->' ,vv ] ;

                else
                    tit=['Results channels ',vv,'->' ,vv1] ;
                end

                hFig = figure(df);
                subplot(1,2,iter);
                title(tit,'Color', 'r')
                set(hFig, 'Position', [10 10 10 10 ]);
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
                text(HRC2+0.01,0.12,['ci_u', ' -- ' , num2str(round(HRC2,1)),' ms']);
                text(HRC2+0.01,0.15,['Delay', ' -- ' , num2str(round(meanHRC,1)),' ms'],'Color','red','FontSize',12);
                text(HRC2+0.01,0.18,['ci_l', ' -- ' , num2str(round(HRC,1)),' ms']);

                grid('on');
                xlabel('Delay');
                ylabel('NRMSE');
                hold off
                
                               
                
        end
                
    end
end
end