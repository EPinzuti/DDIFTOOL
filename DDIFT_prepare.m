
function Data_analysis = DDIFT_prepare(varargin)

% DDIFT_prepare check the configuration input cfg and data structure.
% It prepares the data for subsequent analysis.
% You can call this function as follows: 
%
%  data_prepare=DDIFT_prepare(cfg, data)
%
%
%The function DDIFT_prepare was inspired, and it follows many aspects of
% cheking input parameters, by TEpreapare.m from TRENTOOL toolbox:
%
% "M. Lindner, R. Vicente, V. Priesemann, M. Wibral (2011). TRENTOOL: 
% a Matlab open source toolbox to analyse information flow in time series
% data with transfer entropy. BMC Neurosci., 12, p.119"
%
%
%*** INPUT PARAMETERS  as Fieldtrip format. Row data structure, data must
%                      contains:
%
%     
%       .trial     = cell array (nr of channels x nr of samples) containing
%                    the data for each trial
%       .time      = cell (1xnr of samples) containing the time indices for
%                    each trial (in seconds)
%       .label     = cell (1xnr of channels), containing the labels
%                    (strings) of channels included in the data
%       .fsample   = value of sampling rate (in Hertz)
%
%*** cfg PARAMETER STRUCTURE must contains:
%
%
%      .delay = [start  end  step size]  e.g. [0  60  5]. 
%                 It defines the delay grid of REG.  
%                 The spacing should agree with the system time scale.
%       
%
%      .order = Model order and model type.
%               Integer numbers are required to work with the parametric 
%               model (1,2 ,3,...),'inf' for the non-parametric model
%
%
%      .es  = [dimension , delay].Specify dimension and delay of embedding. 
%             If ‘compute’ cfg.param and cfg.tau have to be 
%             provided by the user
%
%
%      .subsampling = Value to downsample the time series. 
%                     If ‘compute’  the user can provide the input  
%                     from the console once mutual information is performed.
%
%
%      .channel = cell {1 x n.chan} Names of  channels  to analyse(string). 
%
%
%      .combination = manual’ to use a specific channel pairs in the analyze 
%                      or ‘All’ for all possible channels combination.
%                     If cfg.combination is set to ‘manual’, cfg.channel
%                     have to be a  n. channels x 2  cell array.
%
%       .time     = First and last point of the time vector of interest 
%                   (in seconds)  for analysis.
%
%       .sampling  = Sampling rate of time series  1/fsample
%
%*** cfg OPTIONAL
%
%
%       .param     = [min max] Range of embedding dimensions to test 
%
%       .tau       = [min max] Range of embedding delays  to test
%
%       .padding = If ‘True’ a zero-padding  of length d-1 is add 
%                       at the beginning of each trials.
%
%       .numvalidate = Sample points taken from the end  of the time series 
%                      or validation.The value needs to be higher or equal
%                      than the  dimension of embedding.
%
%
%
%       .testing    = True’ to test pairs channel/time series 
%                     in both directions or ‘False’  test only direction as
%                     entered in  ‘channel’ by the user
%
%
%       .display      = ‘True’ or ‘False’   for plotting results
%
%       .verbosity    = Output of console
%
%       .parallel      Set to ‘True’ for parallel computing. It takes 
%                      the default number of workers. 
%
%*** OUTPUT
%
% DATA    = The output of this function is the data from the input with
%           the added structure D^2IFT_prepare.  


%   DDIFT_prepare
%         
%         .pre_ch=structure with selected channels from the user. Specified  
%                 channels are kept for subsequent analysis
%      
%         .timeindices=Indices in samples of the time of interest 
%                      selected by the user [start end ]
%
%         .es=Estimated d and tau by embedding routine
           % is a cell n.channel pair x 4
           %where:
           %the first two columns are names of analysed channel
           %third column direction chan104 --> chan123
           %fourth column direction chan123 --> chan104
           
             %                            A--->B    B-->A
             % { ['chan104'] ['chan123'] ['20 1'] ['18 1']}
             % { ['.......'] ['.......'] ['.. .'] ['.. .']}
           
             % both directions are tested. In principle the direction 
             % with higher d dimension is the most interesting, indeed
             %if the slave carries information of the driver, higher d dimension is
             %required to unfold the geometry of the attractor (skew
             %product embedding)
             
             
%
%         .subsampling=Estimated subsampling values by mutual information
          % is a cell number of unique channel x2 
          % { ['104'] [12]}
          % { ['103'] [11]}
          % { ['..']  [..]}
%         auto-mutual inf is computed once even if the user provides the
%         same channel name multiple times.


%%
if isfield(varargin{1},'es') && isstruct(varargin{1}) && isstruct(varargin{2}) && isfield(varargin{2},'trial')
    cfg =  varargin{1};
    data = varargin{2};
else
    error('\DDIFTOOL: incorrect input values, see help!');
end
%%


LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;
LOG_DEBUG_COARSE = 3;
LOG_DEBUG_FINE = 4;

if ~isfield(cfg, 'verbosity'), cfg.verbosity = 'info_m'; end;
%%

msg = 'Checking data and config';

console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);


 
 % check the data structure
if ~isfield(data, 'trial'),
    fprintf('\n')
    error('DDIFTOOL error: data must be in ''.trial''-field, see help!');
end;
if ~isfield(data, 'time'),
    fprintf('\n')
    error('DDIFTOOL error: data contains no ''.time''-field, see help!');
end;
if ~isfield(data, 'label'),
    fprintf('\n')
    error('DDIFTOOL error: data contains no ''.label''-field, see help!');
end;
if ~isfield(data, 'fsample'),
    fprintf('\n')
    error('DDIFTOOL error: data contains no ''.fsample''-field, see help!');
end;
 %check data using checkdata from Fieldtrip
[data] = ft_checkdata(data, 'datatype','raw');
 % check whether time axes and trials have the same number of entries. Taken
 %from TRENTOOL TE_prepare.m
if iscell(data.time) % one time axis per trial
    for tt=1:size(data.trial,2) % for each trial
        if ~( size(data.time{tt},2) == size(data.trial{tt},2) )
            errorstr=strcat('DDIFTOOL error! incorrect number of samples in time axis or trial detectedin trial Nr:',...
                num2str(tt),...
                ', samples: ',num2str(size(data.trial{tt},2)),...
                ', timeindices: ',num2str(size(data.time{tt},2)) );
            fprintf('\n')
            error(errorstr)
        end
    end
else % time is a single vector
    for tt=size(data.trial,2)    % for each trial
        if ~( length(data.time) == size(data.trial{tt},2) )
            disp('in trial Nr: ')
            disp(num2str(tt))
            mes=strcat('DDIFTOOL error! incorrect number of samples in time axis or trial, detected in trial Nr:',num2str(tt));
            error(mes);
        end
    end
end

%check configuration and set defaults

% set sampling rate from data
 %data.fsample;
if ~isfield(cfg, 'subsampling') && strcmp(cfg.subsampling,'compute')
     msg = 'specify a value or ''compute'' if a subsampling migth be necessary (reccomanded)';

     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
     if strcmp(cfg.subsampling,'compute')
        cfg.subsampling =[];
     end
     
    
     
end;


if ~isfield(cfg, 'numvalidate')
     cfg.numvalidate = [];    
     msg = 'default parameter for numvalidate will be 1/5 of the time series length';

     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
     
end;
if ~isfield(cfg, 'testing')
     cfg.testing = 'True';    
     msg = 'default parameter for testing=True';

     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
     
end;
if ~isfield(cfg, 'display')
     cfg.display = 'True';    
     msg = 'default parameter for display=True';

     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
     
end;

if ~isfield(cfg, 'delay')
    fprintf('\n')
    error('you need to specify delay candidates');
else
    start=cfg.delay(1);
    end_d=cfg.delay(2);
    step_s=cfg.delay(3);
    cfg.delay=linspace(start,end_d,(end_d/step_s)+1);
    
end
if ~isfield(cfg, 'time')
    fprintf('\n')
    error('you need to specify time of interest');
end   

if ~isfield(cfg, 'order')
    fprintf('\n')
    error('you need to specify model order');
end
if ~isfield(cfg, 'es')
    fprintf('\n')
    error('you need to specify embedding');
end
if ~isfield(cfg, 'channel')
    fprintf('\n')
    error('you need to specify which channels to analyse');
end
if ~isfield(cfg, 'combination')
    fprintf('\n')
    error('you need to specify how to combine channels for testing');
end

if ~iscell(cfg.channel)
    
    fprintf('\n')
    error('channel needs to be a cell ');
end
if ~isfield(cfg, 'padding')
     cfg.padding = 'False';    
     msg = 'default parameter for Zero-padding=False';

     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
end
if ~isfield(cfg, 'parallel')
     cfg.parallel = 'False';    
     msg = 'default parameter for parallel computing=False';

     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
end
if ~isfield(cfg, 'fold_name')
       
     error('Please provide a foleder name to save data analysis ');
end

%user can choose trials to analysise 
% if ~isfield(cfg, 'trial') 
%         
%      msg = 'default all trials';
% 
%      console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
%      trial_set=0;
% elseif isempty(cfg.trial)
%     msg = 'default all trials';
% 
%     console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
%     trial_set=0;
% 
%     
% end


%%


% check if channel or channelcombinations are defined

msg = 'constructing data structure for analysis of different channels';
console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
if size(cfg.channel,2)==2
    
    compare=reshape(cfg.channel,1,size(cfg.channel,1)*size(cfg.channel,2));
    number_channel=length(compare);
else
    
    
    if size(cfg.channel,1)>size(cfg.channel,1)
       compare=unique(cfg.channel);
       number_channel=length(compare);
    else
        cfg.channel=cfg.channel';
        compare=unique(cfg.channel);
        number_channel=length(compare);
      
    end
end

%cfg.channel need to be a cell list

%checking if same name is provided as in data struct
    
for ch=1:number_channel
   
    

    if ismember(compare{ch},data.label)==0
        
        
        fprintf('\n')
        error('DDIFTOOL error ! channel name does not match with channel name in data');
    else
        continue
    end
end
%take time of interest
if size(cfg.time,1) > 2 || size(cfg.time,2) >2
    fprintf('\n')
    error('DDIFTOOL error! cfg.time has more than two entries');
end
if size(cfg.time,1)>size(cfg.time,2)
    cfg.time=cfg.time';
end
% read time values of the data
if iscell(data.time)
    alltime=cell2mat(data.time(1));
else
    alltime=data.time;
end

% find correct indices for the samples 
% to be used later

timeindices=zeros(1,2);

for ii = 1:size(cfg.time,2)
    [col]=nearest(alltime, cfg.time(ii));
    timeindices(ii)=col;
end


% check indices 
if timeindices(1) >= timeindices(2)
    error(['DDIFTOOL error ! Something seems to be wrong with your ' ...
        'time indices %d and %d: time index 1 >= time index 2!'], ...
        timeindices(1), timeindices(2))
else
    DDIFT_prepare.timeindices = timeindices;
   
end     

% time_s=cfg.time(1);
% time_e=cfg.time(2);
%if all channel need to be tested 
if strcmp(cfg.combination,'All') 
    
% if trials are not of the same lengths gives error. XXX
%it is better to define again which trials to analyse using fieldtrip.


%Subset of channels have to be tested with all combination.

    [idx_ch,~,~]=channel_select(data,cfg);
    number_trial=length(data.trial);
    idx_ch_t=unique(idx_ch);
    data_p=struct();
%     channel=zeros(number_trial,length(data.trial{1,1}));
    for tr=1:number_trial

        temp=data.trial{1,tr};
        for ch=idx_ch_t'
           try
                data_p.channel(ch).channel(tr,:)=temp(ch,:);
            catch ME
               if (strcmp(ME.identifier,'MATLAB:subsassigndimmismatch'))

                    error(['DDIFTOOL error !  ' ...
           ' trials might have different length. Please provide trials ' ...
           'with the same length (you can do it with FieldTrip functions) ']);
               end
                
            end
        end
    end

    DDIFT_prepare.pre_ch=data_p.channel;  
elseif strcmp(cfg.combination,'manual') 
    compose=size(cfg.channel,2);
    if  compose>2
        
%         fprintf('\n')
        
         error('DDIFTOOL error ! you need to provide a list of pairs n x 2')
   
    else
        
        [idx_ch,~,~]=channel_select(data,cfg);
        idx_ch_t=unique(idx_ch);
        number_trial=length(data.trial);
        data_p=struct();
%         channel=zeros(number_trial,length(data.trial{1,1}));
        for tr=1:number_trial

            temp=data.trial{1,tr};
            for ch=idx_ch_t'
                
                try
                    data_p.channel(ch).channel(tr,:)=temp(ch,:);
                catch ME
                    
                    if (strcmp(ME.identifier,'MATLAB:subsassigndimmismatch'))

                        error(['DDIFTOOL error ! ' ...
               ' trials might have different length. Please provide trials ' ...
               'with the same length (you can do it with FieldTrip functions) ']);
                    end
                
                end
                

            end
        end

    %     data_p1= {data_p.channel};
    %     dd=~cellfun(@isempty,data_p1{1,:});

        DDIFT_prepare.pre_ch=data_p.channel; 
    end
else
    fprintf('\n')
        error('DDIFTOOL error ! cfg.combination support "All" or "manual" ');
end

   
% build data structure for analysis with only data points of interest

chg=zeros(length(idx_ch'),1);
t=1;
for jj=idx_ch'
    if isempty(find(jj==chg,1))
    
        DDIFT_prepare.pre_ch(jj).channel= DDIFT_prepare.pre_ch(jj).channel(:,DDIFT_prepare.timeindices(1):DDIFT_prepare.timeindices(2));
        length_series=length(DDIFT_prepare.pre_ch(jj).channel);
        chg(t,:)=jj;
        t=t+1;
    else
        continue
    end
end
% % from now consider only trials specify by the users, if cfg.trial is empty
% % all trials are considered
% 
% if trial_set==0
%     %do nothing all trials are kept
% else
%     
%     % function trial select  XXXX
% end
%     
%% Set Parallel Computing
if strcmp(cfg.parallel,'True') 
 
   % check if it possible to set parallel computing
   [parallel_state, poolobj]=set_parallel(cfg);
   inf_paral=poolobj;
   cfg.par_state=inf_paral;
   if parallel_state==0
       
       msg = 'It is not possible to set parallel computing correctly';

       console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
   end
       
   
else
   cfg.par_state=[];
   
end


 
%% set last cfg parameter
%Final preprocessing step
%%
%%estimate embedding parameter
%if es='compute' the embedding estimation module it is called
%if is in the form [20 1] DDIFT_prepare assume that the parameter are known
%and skip this step.

% if size(cfg.es,2)>2
msg = 'Estimating embedding dimension and delay parameters ';
console_output(cfg.verbosity, msg, LOG_INFO_MAJOR)
        
if strcmp(cfg.es,'compute')
    min_d=cfg.param(1);
    max_d=cfg.param(2);
    tau_l=cfg.tau(1);
    
    if length(cfg.tau)>1
        tau_u=cfg.tau(2);
        esRange={round(linspace(min_d,max_d,(max_d+1)-min_d)),round(linspace(tau_l,tau_u,(tau_u+1)-tau_l))};
       
        
    else
        esRange={round(linspace(min_d,max_d,(max_d+1)-min_d)),tau_l};
        
    end

    
    
    DDIFT_prepare.esRange=esRange;
%         for j=1:size(Neco_prepare.pre_ch,2)
    [idx_ch,indices,list]=channel_select(data,cfg); 
%     DDIFT_prepare.indices=indices;
    %if cfg.combination is set to 'all', every possible combination of channel and both direction
    %are tested with range dimension
    % if cfg.combination is set to 'manual' same
    
%     DDIFT_prepare.es=cell(size(indices,1),2);%zeros(size(indices,1),4)
    u=1;
    DDIFT_prepare.es=[list,cell(size(indices,1),2)];
%     rr=zeros(length(esRange),3);
%     rr1=zeros(length(esRange),3);
    for f_c=1:size(indices,1)
       
        f_ch=indices(f_c,1);
        data1=DDIFT_prepare.pre_ch(f_ch).channel;


        t_ch=indices(f_c,2);
        data2=DDIFT_prepare.pre_ch(t_ch).channel;
        driver=(data1);
        slave=(data2);
        
        
%         for ee=1:size(data1,2)
       % for1:
            % use available trials, this might be a problem computationally
            % It also depends on the length of the time-series
            jk=size(driver,1);
            if jk>5
                jk=5;
            else
            end
            [d0, t0, resid]=embed_cv2(slave(1:jk,:),driver(1:jk,:),esRange,cfg.par_state,cfg.verbosity);
            DDIFT_prepare.es{u,3}=num2str([d0,t0]);
            
            % testing other direction
            [d0, t0, resid]=embed_cv2(driver(1:jk,:),slave(1:jk,:),esRange,cfg.par_state,cfg.verbosity);
            DDIFT_prepare.es{u,4}=num2str([d0,t0]);
%             DDIFT_prepare.es{u,3}(1,2)=t0;
%             DDIFT_prepare.es_resid=cell(size(indices,1),2);
%             DDIFT_prepare.es_resid{1}(u,:)=resid;
%             rr(:,u)=d0;
    
            
%             DDIFT_prepare.es{u,4}(1,2)=t0;
%             DDIFT_prepare.es_resid1=resid;
%             rr1(:,u)=d0;
%         end
        %end
%             Neco_prepare.es(f_c,6)=resid;

        u=u+1;

%   
%                       
    end


    display( DDIFT_prepare.es)
    % solution for multiple input split taken from  Walter Roberson  
    valstring = input('Please specify integer for d and tau:', 's');
    valparts = regexp(valstring, '[ ,]', 'split');
    values_param = str2double(valparts);
    
    if length(values_param)<2 || length(values_param)>2
        
        error('DDIFTOOL error ! you need to provide two inputs') 
    else
    end
    %if user provide a letter regexp gives a Nan
    % check if it contains a zero
    if sum(values_param==0)>0
        zer=1;
    else 
        zer=0;
    end
    
    nn=isnan(values_param);
    rest=fix( values_param) ~=  values_param;
    
     while any(nn)|| any(rest) || zer
         valstring = input('Your input is not an integer or it is a zero. Provide an integer number for d and tau(parameters can not be zero :', 's');
         valparts = regexp(valstring, '[ ,]', 'split');
         values_param = str2double(valparts);
         nn=isnan(values_param);
         rest=fix( values_param) ~=  values_param;
         if sum(values_param==0)>0
            zer=1;
         else
             zer=0;
         end
    end
    
%          
     
    dimens= num2str(values_param(1));
    delay_e= num2str(values_param(2));
    msg = ['Parameters for dimension:  ',dimens, ' and delay: ',delay_e];

    console_output(cfg.verbosity, msg, LOG_INFO_MAJOR)
    %from now cfg.es is set to the user input after estimated parameter
    cfg.es=[values_param(1) values_param(2)];
    %call embedding  module
    
else
    dimens= num2str(cfg.es(1));
    delay_e= num2str(cfg.es(2));
    msg = ['Parameters for dimension:  ',dimens, ' and delay: ',delay_e];

    console_output(cfg.verbosity, msg, LOG_INFO_MAJOR)
    
end

%check if the time-series need to be downsample using auto-mutual
%information

%mutual information is computed with classicla k-neighboor
%k=number of nearest points
%n window of mutual information
%y=time_series
if strcmp(cfg.subsampling,'compute')

     msg = 'Checking if the time series need to be downsampled';


    console_output(cfg.verbosity, msg, LOG_INFO_MAJOR)
    % set path for TSTOOL functions, taken from TRENTOOL (TEarch)
    arch;
    [idx_ch,indices,~]=channel_select(data,cfg); 
    uniq=unique(indices)';
    %DDIFT_prepare.sunsambling is a cell n.unique channel x2 
              % { ['104'] [12]}
              % { ['103'] [11]}
              % { ['..']  [..]}
              
    DDIFT_prepare.subsampling=cell(size(uniq,2),2);
   
    i=1;
%     chg=zeros(length(idx_ch'),1);
    t=1;
    % idx_ch:t carries unique channel, so mutual information is computed
    % once even if the user specified the same channel multiple time
    for kk=idx_ch_t'
           msg = [' computing mutual information '];
            
           console_output(cfg.verbosity, msg, LOG_INFO_MAJOR)
%         if isempty(find(kk==chg))
            DDIFT_prepare.subsampling{i,1}=data.label{kk};
            mi_c=zeros(50,1);
            % set to first trial, it can compute mutual information for
            % each trial (uncomment number_trial).
            %usually mutual information is similar to every trial
            for jk=1%:number_trial
                % use first 6000 if time-series is longer
                if length(DDIFT_prepare.pre_ch(kk).channel(jk,:)) > 6000
                    y=(DDIFT_prepare.pre_ch(kk).channel(jk,1:6000)');
                else
                    y=(DDIFT_prepare.pre_ch(kk).channel(jk,:)');
                end
                    
                       
                [sub,mi]=mutual_subsampling(y,50,3);
                mi_c(:,jk)=mi;
                
%                 DDIFT_prepare.s_samp(i,jk)=sub;
    %       
                
            end
            
            % take average between trials gives a often a smooth function
            % without a clear local minima / not used 
            averag_m=mean(mi_c,2);
            
   
            figure(i)
            %plot average mutual information 
%           plot(  averag_m,'linewidth',2)
%           plot all trials, in this case just one
            plot( mi_c(:,:),'linewidth',1);
            
            ylabel('Information (bits)');
            xlabel(' Time-lag');
           
            title(['Auto-mutual info',data.label{kk}]);
            hold off
            %function taken from David Sampson
%           k=findminima(averag_m); on average 
%        
            k=findminima(mi_c(:,1));
%             hold on
%             chg(t,:)=kk;
            %first entry is first local minimum
            local_mini=num2str(k(1));
            DDIFT_prepare.subsampling{i,2}=k(1);
       
            msg = ['first_minimum mutual info ', data.label{kk},': ' , local_mini];
            
            console_output(cfg.verbosity, msg, LOG_INFO_MAJOR)
%           %save figure automatically in the folder XXX  
%            
            t=t+1;
            i=i+1;
%         else
%             continue
%         end
    end



    prompt = 'specify integer for subsampling.  1 if no subsampling:  ';
    s_v = input(prompt);
  
    while rem(s_v,1) ~= 0
%         msg = 'have to be  an integer:  ';
        prompt1 = 'Your input is not an integer. Specify integer for subsampling.  1 if no subsampling:  ';
        s_v = input(prompt1);
      
%         console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
       
    end    
    % downsampling time-series    
    cfg.subsampling=s_v;
    if s_v>1
        chg=zeros(length(idx_ch'),1);
        i=1;
        for kk=idx_ch'
            %dowsample a channel only one time
            if isempty(find(kk==chg,1))




                DDIFT_prepare.pre_ch(kk).channel=downsample(DDIFT_prepare.pre_ch(kk).channel',s_v)';
                number_points=length(DDIFT_prepare.pre_ch(kk).channel); 
                chg(i,:)=kk;
                i=i+1;
            else
                continue
            end
        end





    else isempty(s_v)
 
        msg = 'No input provided, subsampling is set to 1';

        console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);
        cfg.subsampling=1;

    end
else
    
    if cfg.subsampling>1
         chg=zeros(length(idx_ch'),1);
         i=1;
         for kk=idx_ch'
            
             if isempty(find(kk==chg,1))




                DDIFT_prepare.pre_ch(kk).channel=downsample(DDIFT_prepare.pre_ch(kk).channel',cfg.subsampling)';
                
                number_points=length(DDIFT_prepare.pre_ch(kk).channel); %use of channel to estimate how many points are present after downsampling
                chg(i,:)=kk;
                i=i+1;
             else
                 continue
             end
         end
    else
    end
end


if isempty(cfg.numvalidate)
    
    if cfg.subsampling==1
       n_valid=round(length_series/5);   % put the default value
       cfg.numvalidate=n_valid;
    else
        %if numvalidate is empty and subsampling not one ask user to pass a
        %value
        value_point=num2str(number_points);

       
    
     
         prompt = ['specify integer for numvalidate. After subsamplig there are ', value_point, ' data points available:  '];
         s_v_sub = str2double(input(prompt, 's'));
         while isnan( s_v_sub) || fix( s_v_sub) ~=  s_v_sub
              s_v_sub = str2double(input('Please enter and INTEGER: ', 's'));
         end
        
        
           

        cfg.numvalidate=s_v_sub;
    end
    
else
    cfg.numvalidate=cfg.numvalidate;

end




%%
%add zero padding at the begin of each trials
%make vector of zeros


if strcmp(cfg.padding,'True')   

    padd=zeros(number_trial,cfg.es(1)-1);
    chg_pad=zeros(length(idx_ch'),1);
    i=1;
    for jj=idx_ch'
        % add zero only once to each channel
        if isempty(find(jj==chg_pad,1))
            DDIFT_prepare.pre_ch(jj).channel=[padd,DDIFT_prepare.pre_ch(jj).channel(:,:)];
            chg_pad(i,:)=jj;
            i=i+1;
        else
            
            continue
        end
    end

else
end


%see if all parameter are ok. This mean to check if all the values for validate, embedding and delay are consistent with the number of point available
% to avoid the function to crash
%XXX




%%


DDIFT_prepare.cfg=cfg;

varargin{2}.DDIFT_prepare=DDIFT_prepare;
%make a folder to save data for each subject analysed will be 
%changed to deal with multiple subject
Data_analysis = varargin{2};

dataset_string=cfg.fold_name;

% make a function to save data properly XXX (important for multiple sub
% analyis)
% find where is the core function
pathDDIFT = which('DDIFT_prepare.m');
pathDDIFT = pathDDIFT(1:end-15);
pathfolder = fullfile(pathDDIFT,'data_DDIFT');
full_n1=fullfile(pathDDIFT,'data_DDIFT\');
full_n2=fullfile(pathDDIFT,'data_DDIFT\',dataset_string);
% check if data_DDIFT already exist exist and is a folder
if exist(pathfolder, 'dir')==7
    %just save data_prepared under specific subfolder name es: subject1,2,3
    %check if the user is using the same subfolder name to save data
    %prepared
    if exist(full_n2, 'dir')==7
         m=input('already existing subfolder.Do you want to overwrite it?, y/n (lower case):','s');
       
        
         if strcmp(m,'y') || strcmp(m,'Y')
             
             mkdir(full_n1,dataset_string);
             filename=[full_n2 '/' 'data_prepared'];
             save(filename,'Data_analysis');
%              save(['data_DDIFT/' dataset_string '/data_prepared' ], 'Data_analysis');
         elseif strcmp(m,'n') || strcmp(m,'N')
             m_name=input('Provide a different name (string):','s');
             dataset_string=m_name;
             while exist([full_n1,dataset_string], 'dir')==7
                 
                 m_name=input('Already existing subfolder.Provide a different name (string):','s');
                 dataset_string=m_name;
             end
             %change also conf parameter to subsequent saving
             Data_analysis.DDIFT_prepare.cfg.fold_name=dataset_string;
             mkdir(full_n1,dataset_string);
             full_n2=fullfile(pathDDIFT,'data_DDIFT\',dataset_string);
             filename=[full_n2 '/' 'data_prepared'];
             save(filename,'Data_analysis');

%              save(['data_DDIFT/' dataset_string '/data_prepared' ], 'Data_analysis');
         else
             error('DDIFTOOL error: input not recognized!');
         end
    else
        mkdir(full_n1,dataset_string);
        filename=[full_n2 '/' 'data_prepared'];
        save(filename,'Data_analysis');
%         save(['data_DDIFT/' dataset_string '/data_prepared' ], 'Data_analysis');
    end
else
    %if not create the folder for the first time and first subject
    full_n=fullfile(pathDDIFT,'data_DDIFT\');
    mkdir(full_n,dataset_string);
    filename=[full_n2 '/' 'data_prepared'];
    save(filename,'Data_analysis');
%     save(['data_DDIFT/' dataset_string '/data_prepared' ], 'Data_analysis');
end
msg = 'Data are ready for delay analysis';
console_output(cfg.verbosity, msg, LOG_INFO_MAJOR);

%%




    