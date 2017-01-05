function [path2TSTOOL] = arch
% This function is taken from TRENTOOL (TEarch.m) and modified to add TSTOOL to the 
% path to work with DDIFTOOL
% FUNCTION ARCH
%
% Adds the folder with TSTool functions to the MATLAB search path. The
% functions checks which OS is running and whether architecture is 32 or 64
% bit.
%
% 


% check if TSTool functions are already on the path
% h = which('nn_prepare.m');
% if ~isempty(h)
%     return
% end

path2DDIFT = which('DDIFT_prepare.m');
path2DDIFT = path2DDIFT(1:end-15);
path2TSTOOL = fullfile(path2DDIFT,'TSTOOL_functions');

extension = mexext;		% gets the mex extension expected by the current architecture

switch extension
  case 'mexw64'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_wind', 'mex64');
  case 'mexw32'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_wind', 'mex32');
  case 'mexmaci64'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex mac', 'mex64');
  case 'mexmaci32'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex mac', 'mex32');
    case 'mexa64'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex linux', 'mex64');
  case 'mexglx' 
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex linux', 'mex32');
  otherwise
    error(['D^2IFTOOL ERROR: System not supported! Expected extension: ''' mexext '''']);
end

addpath(path2TSTOOL);       % this folder contains the help files for TSTool mex functions
addpath(path2TSTOOL_mex);   % this folder contains the mex functions themselves
