
function [idx_channel,indices,list_ind]=channel_select(data,cfg)

% output indices fro testing different channels time series
if size(data.label,1)==1
    data.label = data.label';
end

allchannels=data.label;

% specific combination or all channel combination, this is specified in
% cfg.channel

user_ch=cfg.channel;
if strcmp(cfg.combination,'All') 
%no specific combination but specific channel name to be tested
    idx_channel=zeros(length(user_ch),1);
    list={};
    for jj=1:length(user_ch)

        for tt=1:length(allchannels)


            if strcmp(user_ch{jj},allchannels{tt})==1

                idx_channel(jj)=tt;
                list{jj}=allchannels{tt};
            else
                continue
            end
        end


    end
    list=list';
    indices=nchoosek(idx_channel,2);
    list_ind=nchoosek(list,2);
elseif  strcmp(cfg.combination,'manual') 
    
    idx_channel=zeros(size(user_ch,1),2);
    list={};
    for jj=1:size(user_ch,1)
        
        
        for tt=1:length(allchannels)


            if strcmp(user_ch{jj},allchannels{tt})==1

                idx_channel(jj,1)=tt;
                list{jj,1}=allchannels{tt};
            else
%                 continue
            end
            if strcmp(user_ch{jj,2},allchannels{tt})==1

                idx_channel(jj,2)=tt;
                list{jj,2}=allchannels{tt};
            else
%                 continue
            end
            
        end


    end
    indices=idx_channel;
    list_ind=list;
    idx_channel=reshape(idx_channel,size(user_ch,1)*size(user_ch,2),1);
else
end

   


            
    
    