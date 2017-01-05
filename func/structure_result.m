
function [Result]=  structure_result(varargin)


data =  varargin{1};
cfg1=varargin{2};    

Result=struct('dlResult',[],'est_param',[],'cfg',[]);
Result.dlResult.channel.test.result=struct('NRMSE',{});
Result.dlResult.channel.test.result=zeros(length(data.DDIFT_prepare.cfg.delay),6);
Result.est_param.channel.test=struct();
Result.cfg=cfg1;
