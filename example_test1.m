%example test
% Below, two examples of the toolbox. The first one reproduce the
% LORENZ_ROSSLER figure 4 of the paper "A statistical framework to infer
% delay and direction of information flow from measurements of complex system
%"(Schumacher et al.,2015).
% The second, is an example of analysis with neural data of a monkey
% (reproduce figure 4.6 of the thesis) to show pre-procesing and other fuctionalities
% implemented in the toolbox.
%
%INSTALLATION
% 
% Just add the folder to the matlab path. Use Make file to compile. If you
% compile cholesky decomposition is performed with a c function coded by
% Rasmussen, otherwise the toolbox use a matlab version automatically.
% TSTOOL functions are added to the path automatically so you do not need to
% addd them.

%COMMENT 
% I tested the toolbox with matlab 2016a on windows and matlab 2014a on
% linux. 
%With cfg.parallel='True' the number of workers used is the default
%(#function plot_connectivity.m requires adjustments#, but it works.)
%XXX
% check if there are enough data points for delay candidates and parameters
% chooesed, otherwise statistical model crash ( need to be done in DDIFT_prepare.m)


%%
%%  Lorenz_rossler system 
% see where is the core function and data
clear all
pathDDIFT = which('DDIFT_prepare.m');
pathDDIFT = pathDDIFT(1:end-15);
%simply load data
str=strcat(pathDDIFT,'dlLorRoeWienerdt0.mat');
load(str);

fake_data=reshape(dlLorRoeWienerdt0,6,9000)';
driverM1=sum(fake_data(1001:6000,1:3),2);
slaveM1=sum(fake_data(1001:6000,4:6),2);

%prepare data in Fieldtrip format

datt=[driverM1';slaveM1'];
time=linspace(0,1,5000);
data=struct();
data.fsample=0.001;
data.label={'A1','A2'};%
data.trial{1,1}=datt;
data.time{1,1}=time;

%prepare cfg
cfg=struct();
cfg.delay=[0 60 5]  % provide [start end  step size]    
cfg.order=3;
cfg.verbosity='info_m';
cfg.es=[20 1] 
cfg.param=[5 22];%test range of dimension and tau only if cfg.es='compute'
cfg.tau=[1 2];
cfg.subsampling=1%'compute';    %  mutual information using opentstool if 'compute'  see later
cfg.numvalidate= 1000;   
cfg.display='True'; %make a plot
cfg.testing='False'; % if 'true' testing both causal and acausal hypothesis
cfg.channel={'A1','A2'};%
cfg.combination='manual'; %if 'All' ceck for all possible pairs combination if 'manual' the user privide specific pairs of channels    % XXX correct for manual done
cfg.time=[0 5] % time of interest (to define a specific window)
cfg.sampling=0.0001
cfg.parallel='True' %use parallel computing
cfg.padding= 'False'% zero-padding
cfg.fold_name='EX_Lorenz'

% prepare data
data_prepare=DDIFT_prepare(cfg,data);

[output]=DDIFT_delays_analysis(data_prepare);

REG_plot(data_prepare,output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Monkey data for an example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
% unzip the folder and do preprocessing or load data immediately,below.
% filen=which('visual_grating_monkey_b.zip');
% unzip(filen);

% pathDDIFT = which('TimeStamp-ECoG-100803-2.mat');
% %pathDDIFT = which('TimeStamp-ECoG-101019-1.mat');
% pathDDIFT1 = pathDDIFT(1:end-27);
% %load time
% Y=load(pathDDIFT);
% %simply load data
% str=strcat(pathDDIFT1,'ECoG-2.mat');
% load(str)
% %load lay
% str1=strcat(pathDDIFT1,'lay.mat');
% load(str1)
% 
% %preprocessing
% data.trial=X(1:128,:);
% data.time={Y.X(1,:)};
% data.fsample=1000;
% data.label=lay.label(1:128);
% data.trial={data.trial};
% data.trial=double(data.trial{1});
% data.trial={data.trial};
% cfg=[];
% cfg.layout=lay;
% layout = ft_prepare_layout(cfg, data_1);
% %plot layout
% cfg.layout = layout; 
% ft_layoutplot(cfg,[],1);

% %EventData;
% trigger =X(129,:); 
% sample  = Y.X;
% pretrig  = -data.fsample*2;
% posttrig =  data.fsample*2;
% 
% %data each orientation in the experiment based on trigger
% orientation=45;
% trl = [];
% for j = 2:(length(trigger)-2)
%   trg1 = trigger(j);
%   trg2 = trigger(j+1);
%   trg3 = trigger(j+1);
%   trg4 = trigger(j-1);
%   
%   if orientation==45
%       if trg1 > 400 && trg2 < 750 && trg3 >= 650 && trg4 <= 400      
%         trlbegin = sample(j) + pretrig;       
%         trlend   = sample(j) + posttrig;       
%         offset   = pretrig;
%         newtrl   = [trlbegin trlend offset];
%         trl      = [trl; newtrl];
%       end
%   elseif orientation==90
%       if trg1 > 400 && trg2 < 1050 && trg3 >= 950 && trg4 <= 400      
%         trlbegin = sample(j) + pretrig;       
%         trlend   = sample(j) + posttrig;       
%         offset   = pretrig;
%         newtrl   = [trlbegin trlend offset];
%         trl      = [trl; newtrl];
%       end
%   elseif orientation==135
%       if trg1 > 400 && trg2 < 1400 && trg3 >= 1300 && trg4 <= 400      
%         trlbegin = sample(j) + pretrig;       
%         trlend   = sample(j) + posttrig;       
%         offset   = pretrig;
%         newtrl   = [trlbegin trlend offset];
%         trl      = [trl; newtrl];
%       end
%   else
%   end
%      
% end
% 
% %define trial
% cfg=[];
% cfg.trl=trl;
% data_1=ft_redefinetrial(cfg,data);
% save('monkey_data','data_1')
clear all
load('monkey_data')




cfg=struct();
cfg.delay=[0 10 1]  % provide [start end  step size]    of delay Delay grid for REG. Spacing should agree with system
                                                         %time scales to be meaningful
cfg.order=2; % volterra order
cfg.verbosity='info_m';
cfg.es=[25 1] %'compute'% %if set to 'compute' provide param and tau to test d and tau
cfg.param=[22 26];   
cfg.tau=[1 ];
cfg.subsampling=11;%'compute'; %   / compute mutual inf with tstool functions (it uses only one trial but it can be computed for each trial)
cfg.numvalidate=40 ;   % for validation if empty you will provide it from the console
cfg.display='True'; %make a plot
cfg.testing='False'; % if 'true' testing both causal and acausal hypothesis
cfg.channel={'chan123','chan104'};%
cfg.combination='All'; %if 'All' check for all possible pairs combination if 'manual' the user provides specific pairs of channels    % XXX correct for manual done
cfg.time=[0 2];% time window to analyse
cfg.sampling=0.001;
cfg.parallel='True'; %use parallel computing
cfg.padding= 'True';% zero-padding
cfg.fold_name='EX_Monkey'; % where to save data

data_prepare=DDIFT_prepare(cfg,data_1);
[output]=DDIFT_delays_analysis(data_prepare);

output.cfg.layout=lay;
output.cfg.export=1;

plot_connectivity(output, data_prepare)
REG_plot(data_prepare,output)
