
function plot_connectivity(result,data)

%This function is inspired and takes some plotting routines (head layout) 
%from TEplot2D.m (TRENTOOL toolbox).

% It allows to plot the head layout (see layout Fieldtrip) and arrows between
% analysed channel that showed an interaction, color coded as inverse NRMSE
% ( percentage of reconstruction). 
% 
%
%    -directed information transfer
%    -bidirectional information transfer
%    -putative synchronization (see Schumacher et al.,2015) is not shown.
%
%%%INPUT 
%
%   data= the initial structure in the Fieldtrip raw format

%   Result= output of DDIFT_delays_analysis.m
%
%   Result contains cfg structure. To work whith the function and change
%   plotting parameter the followiing fields must be added:
%
%   Result.cfg
%
%         .layout=Structure  from Fieldtrip  routine 
%                (see FieldTrip reference for  prepare_layout)
%                 or name of an ascii file *.lay
%
%         .showlabels=False” to plot only analyzed channel labels.  
%                     If “True” all channel labels are displayed.
%
%         .export=1 to export a table with summary of results, 0 no export.   
%
%         .plothead=0 to use FieldTrip function (ft_plot_lay.m) with related layout,  
%                    1 to draw an head layout online.
%
%
% to change plotting parameter:
%
%         
%           .cfg.hcolor        = Color of head cartoon (default = [0,0,0])
%           .cfg.hlinewidth    = Linewidth of the  head plot(default = 2)
%           .cfg.emarker       = Marker symbol (default = 'o')
%           .cfg.ecolor        = Marker color (default = [0 0 0] (black))
%           .cfg.emarkersize   = Marker size (default = 2)
%           .cfg.efontsize       = Font size of electrode labels/numbers (default = 8 pt)
%                                 when cfg.electrodes = 'numbers' or 'labels'
%           .cfg.efontcolor      = Font color of electrode labels/numbers when
%                                 cfg.electrodes = 'numbers' or 'labels'
%                                 (default = [0 0 0])
%           .cfg.hlmarker        = Highlight marker(default = 'o')
%           .cfg.hlcolor         = Highlight marker color (default = [1 0 0] (red))
%           .cfg.hlmarkersize    = Highlight marker size (default = 4)
%     
%
%Color of arrow based on inverse nmrse-
% For each iteration draw an arrow if there is an estimate delay
%%%%%
%15/11/16 If testing='False' in the non-tested direction there is not a 0
% this gave an error. Corrected



%% plot sensor level
%check the data structure and check default parameter

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

if ~isfield(result, 'dlResult'),
    fprintf('\n')
    error('DDIFTOOL error: result struc from DDIFT_delay_analysis');
end;
if ~isfield(result, 'est_param'),
    fprintf('\n')
    error('DDIFTOOL error: result struc from DDIFT_delay_analysis');
end;

if ~isfield(result.cfg,'channelcombination')
    error('DDIFTOOL error: Please provide list of channels used in the causal analysis".');
end;
if ~isfield(result.cfg,'layout')
    error('DDIFTOOL error: Please provide layout from fieltrip".');
end;

if ~isfield(result.cfg,'test')
    error('DDIFTOOL error: Please provide test".');
end;
if ~isfield(result.cfg,'iteration')
    error('DDIFTOOL error: Please provide testing".');
end;
if ~isfield(result.cfg,'delay')
    error('DDIFTOOL error: Please provide delay ".');
end;

if ~isfield(result.cfg,'export');     result.cfg.export = 0;      end;

if ~isfield(result.cfg,'showlabels'); 
    result.cfg.showlabels = 'false';

end;
cfg=result.cfg;
% set default
%Taken form from TEplot2D.m %(TRENTOOL toolbox).
if ~isfield(cfg,'emarker');      cfg.emarker = 'o';     end;
if ~isfield(cfg,'ecolor');       cfg.ecolor = [0 0 0];  end;
if ~isfield(cfg,'emarkersize');  cfg.emarkersize = 2;   end;
if ~isfield(cfg,'efontsize');    cfg.efontsize = get(0,'DefaultAxesFontSize');end;
if ~isfield(cfg,'efontcolor');   cfg.efontcolor = [0 0 0];  end;
if ~isfield(cfg,'arrowcolorneg');   cfg.arrowcolorneg  = [0 0 1];     end;
if ~isfield(cfg,'plothead');     cfg.plothead = 0;      end;
if ~isfield(cfg,'hlmarker');     cfg.hlmarker = 'o';    end;
if ~isfield(cfg,'hlcolor');      cfg.hlcolor = [1 0 0]; end;
if ~isfield(cfg,'hlmarkersize'); cfg.hlmarkersize = 4;  end;

if ~isfield(cfg, 'hcolor');        cfg.hcolor = [0 0 0];     end;
if ~isfield(cfg, 'hlinewidth');    cfg.hlinewidth = 2;       end;

%%
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;
msg = 'Preparing connectivity plot';

console_output(data.DDIFT_prepare.cfg.verbosity, msg, LOG_INFO_MAJOR);

% check whether the entry in cfg.layout is a structure (if not assume ascii
% file)
%Taken form from TEplot2D.m %(TRENTOOL toolbox).
if isstruct(cfg.layout) && all(isfield(cfg.layout, {'pos';'label'}))
    lay = cfg.layout;
    xcoord = lay.pos(:,1);
    ycoord = lay.pos(:,2);
    allchannels = lay.label;
    index = 1:size(xcoord,1);

else
    % read .lay file
    filename = strcat(cfg.layout,'.lay');
    [index,xcoord,ycoord,zcoord,v5,allchannels] = textread(filename,'%d %f %f %f %f %s',-1);
end;

% Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45 this is
% needed to make a head on fly without calling ft_plotlayout
%Taken form from TEplot2D.m %(TRENTOOL toolbox).
if cfg.plothead==1 
    xcoord = 0.8*((xcoord-min(xcoord))/(max(xcoord)-min(xcoord))-0.5);
    ycoord = 0.8*((ycoord-min(ycoord))/(max(ycoord)-min(ycoord))-0.5);
end
% matrix with used channel in the causal analysis 
channel_com=cfg.channelcombination;

%find position of used channel from layout 
indices_com=zeros(size(channel_com,1),2);
for ii=1:size(channel_com,1)
    [ch, idxx]=ismember(channel_com(ii,1),lay.label);
    [chz, idzz]=ismember(channel_com(ii,2),lay.label);
    indices_com(ii,1)=[idxx];
    indices_com(ii,2)=[idzz];
           
end
ch_used_name=unique(channel_com);

%use the indices to find the position coordinates
chan_uniq=unique(indices_com)';
chan_xpos=zeros(size(chan_uniq,1),1);
chan_ypos=zeros(size(chan_uniq,1),1);
for jj=1:length(chan_uniq)
    
    chan_xpos(jj,:)=xcoord(chan_uniq(jj),1);
    chan_ypos(jj,:)=ycoord(chan_uniq(jj),1);
end

% find which pair channel have to be be connected  with an arrow
all_result=[];
all_result.channel=cell(1,1);
for i=1:cfg.iteration
    for k=1:cfg.test
            all_result.channel{i}(1,k)=result.est_param.channel(i).test(k).param_final_delay;
            if cfg.test==1
                all_result.channel{i}(1,2)=0;
            end
                
    end
end


table_r=channel_com;
c_del=cell(size(table_r,1),3);
interaction=cell(size(table_r,1),2);
for ii=1:size(table_r,1)
    check_d=all_result.channel{ii};  % check if for the channel test pair there are only zero (mean no interaction )
    if ~any(check_d)
        % Putative synchronization is not used anymore because could give
        % wrong result.
        % % fitting a line to the result, later ask user in case of ambigous estimation
        
%         if cfg.test==2
%             
%         for ij=1:2    
%             len_ff=length(result.est_param.channel(ii).test(ij).param_pNrmse_t);
%             x=linspace(0,cfg.delay,len_ff)';
%             y=result.est_param.channel(ii).test(ij).param_pNrmse_t;
%             [p S]=polyfit(x,y,1);
%             y1=polyval(p,x);
% %             figure(ii)
% %             subplot(1,2,ij)
% %             plot(x,y,'o',x,y1)
% %             ylim([0 1])
%             if S.normr>4
%                 interaction{ii,ij}='putative sync'; %XXX
%                 %both direction should not differ much
%             else
%                 interaction{ii,ij}='no interaction';
%             end
%             %double check both direction should matchXXX
%             
%         end
%         else
%         end
        ddelay=0;
        
        c_del{ii,1}=table_r{ii,1};
        c_del{ii,2}=table_r{ii,2};
        c_del{ii,3}=ddelay;
    elseif (any(check_d==0))~=1 
         % check if there is at least one zero, if not bidirectional (both
         % direction show a delay values)
         
         %XXX finish
         
         interaction{ii,1}='bidirectional'; % possible scenario in case of weak coupling
         interaction{ii,2}='bidirectional';
    else
        

        pos_del=find(all_result.channel{ii}~=0);
        ddelay=all_result.channel{ii}(1,pos_del);
        c_del{ii,3}=ddelay;
        if pos_del==2
% 
            c_del{ii,1}=table_r{ii,1};
            c_del{ii,2}=table_r{ii,2};
            interaction{ii,2}='directed inf flow';
            interaction{ii,1}='no_interaction';
        else
            c_del{ii,1}=table_r{ii,1};
            c_del{ii,2}=table_r{ii,2};
            interaction{ii,1}='directed inf flow';
            interaction{ii,2}='no interaction';
        end
       
    end
   
end

new_table=[c_del,interaction];  % add confidence interval and make this as real table that can be exported

% the table is not ready yet there could be interactions but not discernible
% delay can be found from the reconstruction, like sinchronization
% if not all value are 1 the can be checked trying to fit a line to the
% function (not used)


% XXX
legendInfo={}; 

    for i=1:size(new_table,1)

        if new_table{i,3}~=0 || strcmp(new_table{i,4},'putative sync') || strcmp(new_table{i,4},'directed inf flow')
            delay=new_table{i,3};
            label_chan=new_table{i,1};
            label_chan2=new_table{i,2};
            legendInfo{i} = [label_chan '-' label_chan2 ' ' num2str(delay) ' ms' ' ' new_table{i,4}];
        elseif new_table{i,3}~=0 || strcmp(new_table{i,4},'putative sync') || strcmp(new_table{i,5},'directed inf flow')
            delay=new_table{i,3};
            label_chan=new_table{i,2};
            label_chan2=new_table{i,1};
            legendInfo{i} = [label_chan '-' label_chan2 ' ' num2str(delay) ' ms' ' ' new_table{i,5}];
           
        else
            %bididirectional XXX
        end
    end



% Define the outline of the head:
%%
%Taken form from TEplot2D.m %(TRENTOOL toolbox).
figure
if cfg.plothead==1
    % Plot head, ears, and nose:
    rmax=.5;
    l     = 0:2*pi/100:2*pi;
    tip   = rmax*1.15; base = rmax-.004;
    EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
    EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];

    hold on
        
    plot(cos(l).*rmax, sin(l).*rmax, 'color', cfg.hcolor , 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
    
    plot([0.18*rmax;0;-0.18*rmax], [base;tip;base], 'Color', cfg.hcolor , 'LineWidth', cfg.hlinewidth);
    plot( EarX, EarY, 'color', cfg.hcolor , 'LineWidth', cfg.hlinewidth)
    plot(-EarX, EarY, 'color', cfg.hcolor , 'LineWidth', cfg.hlinewidth)
    %hold off
    hold on 
end
%%
if cfg.plothead==0  % use field trip to plot layout in case of no human brain to maintain correct position of electrodes
    
%     ft_layoutplot(cfg, data, 0) % better to use lay plot and put false to label XXX
     if strcmp(cfg.showlabels,'True')
         
         
        ft_plot_lay(lay, 'point', true, 'box', true, 'label', true, 'mask', true, 'outline', true);
     else
         ft_plot_lay(lay, 'point', true, 'box', true, 'label', false, 'mask', true, 'outline', true);
     end
     hold on
     
end
%
for ii = 1:size(xcoord,1)
        plot(xcoord(ii),ycoord(ii),cfg.emarker,'Color',cfg.ecolor,'Linewidth',cfg.emarkersize)
        
end

 
 for jj = 1:size(chan_xpos,1)
        plot(chan_xpos(jj),chan_ypos(jj),cfg.hlmarker,'Color',cfg.hlcolor,'Linewidth',cfg.hlmarkersize)
end

    
    
    
%plot label of used channel    
for ii=1:length(chan_uniq)
        %h=text(channelxcoords(ii)+0.01,channelycoords(ii)+0.01,channelnames(ii));
        h=text(chan_xpos(ii)+0.03,chan_ypos(ii)+0.03,ch_used_name(ii));
        set(h,'FontSize',cfg.efontsize,'Color',cfg.efontcolor);
end

% define type of interaction

con_min = 0; % limit nrmse
con_max = 1;
map = brewermap(100,'YlOrRd');p=map;
% p=colormap('parula')

% p=colorbar; colormap('jet'); %need to be changed
reconstruct=cell(size(new_table,1),1);  % nrmse in percentage
for i=1:size(new_table,1)
    
    if strcmp(new_table(i,4),'directed inf flow')  | strcmp(new_table(i,4),'putative sync') |strcmp(new_table(i,4),'bidirectional')
        
        % table has information of direction / putative sync is always in
        % both column 4 and 5 same as bidirectional
         X_temp=xcoord(indices_com(i,1),1);  
         Y_temp=ycoord(indices_com(i,1),1);
        
         
         X_temp1=xcoord(indices_com(i,2),1);
         Y_temp1=ycoord(indices_com(i,2),1);
        
         X=[X_temp;X_temp1];  % if direction is in pos 1 no swap direction X Y
         Y=[Y_temp;Y_temp1];
         
         if strcmp(new_table(i,4),'putative sync')

           inv=sum(result.dlResult.channel(i).test(1).result(:,2))/length( result.dlResult.channel(1).test(1).result(:,2));
           A=(1-inv);
         elseif strcmp(new_table(i,4),'bidirectional')
             inv=sum(result.dlResult.channel(i).test(1).result(:,2))/length( result.dlResult.channel(1).test(1).result(:,2));
             inv1=sum(result.dlResult.channel(i).test(2).result(:,2))/length( result.dlResult.channel(i).test(2).result(:,2));
             A_t=(1-inv);
             A_t1=(1-inv1);
             A=[A_t;A_t1];
         elseif strcmp(new_table(i,4),'directed inf flow')
             inv=sum(result.dlResult.channel(i).test(1).result(:,2))/length( result.dlResult.channel(i).test(1).result(:,2));
             A=(1-inv);
                 
%          elseif strcmp(new_table(i,5),'directed inf flow')
%              inv=sum(result.dlResult.channel(i).test(2).result(:,2))/length( result.dlResult.channel(i).test(2).result(:,2));
%              A=(1-inv)*100;
            
         else 
             
         end
  
         ind = round((size(p,1)-1)*A(:,1)/(con_max-con_min))+1;

         
         for kk = 1:size(X,2)
            %Taken from TRENTOOL toolbox and modified here
            N_arrow(X(1,kk),X(2,kk),Y(1,kk),Y(2,kk),3,cfg.arrowcolorneg,1,ind,new_table(i,4));
         
            
         end
%          B=num2str(round(A*100));
%         
%          reconstruct{i,:}=[B(1,:) '%' ' - ' B(2,:) '%'];
%         
         if strcmp(new_table(i,4),'bidirectional')
            B=num2str(round(A*100));
        
            reconstruct{i,:}=[B(1,:) '%' ' - ' B(2,:) '%'];
         end
         if strcmp(new_table(i,4),'bidirectional')==0    
            reconstruct{i,:}=[num2str(round(A*100)) '%'];   
         end
    elseif  strcmp(new_table(i,5),'directed inf flow') 
        %if in table column 5 swap channel
         X_temp=xcoord(indices_com(i,1),1);
         Y_temp=ycoord(indices_com(i,1),1);
        
         
         X_temp1=xcoord(indices_com(i,2),1);
         Y_temp1=ycoord(indices_com(i,2),1);
    
         X=[X_temp1;X_temp];  % if direction is in pos 2 swap direction X Y
         Y=[Y_temp1;Y_temp];  

        
         inv=sum(result.dlResult.channel(i).test(2).result(:,2))/length( result.dlResult.channel(i).test(2).result(:,2));
         A=(1-inv);
         
         
         for kk = 1:size(X,2)
            % make arrow to connect channel of interest
            N_arrow(X(1,kk),X(2,kk),Y(1,kk),Y(2,kk),3,cfg.arrowcolorneg,2,ind(kk),new_table(i,5));
            

         end
         reconstruct{i,:}=[num2str(round(A*100)) '%'];
    else 
       reconstruct{i,:}=[]; 
        
    end
    
    
end
legendInfo=legendInfo(~cellfun('isempty',legendInfo)); 
if isempty(legendInfo)==0
    [a, obj]=legend(legendInfo);
    b=findobj(obj,'type','line');
    set(b,'visible','off');
else
end
if cfg.export==1
    export_table=[new_table,reconstruct];
    T= cell2table(export_table,...
        'VariableNames',{'Channel_A' 'Channel_B' 'delay' 'Type_Interaction_A_B' 'Type_Interaction_B_A' 'reconstructibility'});
%     save('table_result','T');
    %open directory core function and subfolder data_DDIFT
    pathDDIFT = which('DDIFT_prepare.m');
    pathDDIFT = pathDDIFT(1:end-15);
    full_n2=fullfile(pathDDIFT,'data_DDIFT\',data.DDIFT_prepare.cfg.fold_name);
    % or use cd to open the folder
    filename=[full_n2 '\' 'table_of_results'];
    save(filename,'T');
    
    
    
else
end
end
%% plot source level 


