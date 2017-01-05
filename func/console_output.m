
function console_output(verbosity, message, loglevel, varargin)


% print command execution on command lines =message from  (verbosity request)


%verbosity is a mainly a string that describe the execution of the toolboox


% logging levels:
%   'none'         - 0 - no output at all
%   'info_m'       - 1 - major program execution steps (e.g. data 
%                        preprocessing, interaction delay reconstruction, 
%                        TE estimation)
%   'info_sub_rout'- 2 - minor program execution steps (e.g. call to 
%                        subroutines like TEchannelselect)
%   'debug_coarse' - 3 - coarse debug information
%   'debug_fine'   - 4 - very detailed debug information


% get current stack for line info in output 
%dbstack displays the line numbers and file names of the function calls that led to the current pause condition, listed in the order in which they execute
stack = dbstack;
if length(stack) > 1
    stack = stack(2:end);
else % if function is called directly
    stack.file = 'base';
    stack.name = 'base';
    stack.line = nan;
end




switch verbosity     % verbosity is assignet to 0- to 4 and with looglevel is checked if this level requested is the same as loglivel if yes then is printed
    case 'none'
        verbosity = 0;
    case 'info_m'
        verbosity = 1;
    case 'info_sub_rout'
        verbosity = 2;
    case 'debug_coarse'
        verbosity = 3;
    case 'debug_fine'
        verbosity = 4;
    otherwise
        fprintf('\n')
        error(' Unknown verbosity level!')
end

if loglevel > verbosity 
    return
end

if iscell(message)
    message = cell2txt(message, length(stack));
    
    if ~isempty(varargin)
        message = [varargin{1} message];
    end
end


if loglevel == 1
    fprintf(['\n\n%s - line %d: \n' message '\n'], stack(1).file, stack(1).line);     
        
elseif loglevel == 2 || loglevel == 3
    fprintf('\n')
    for i=1:length(stack)-1;
        fprintf('   ')
    end
    fprintf(['%s - line %d: ' message], stack(1).file, stack(1).line);
    
elseif loglevel == 4
    for i=1:length(stack)
        fprintf('%s - line %d\n', stack(i).file, stack(i).line);
        for j=1:i
            fprintf('\t')
        end
    end
    fprintf(['msg: ' message '\n']);
    
end
end


function txt = cell2txt(msg, indent)

[rows, cols] = size(msg);
txt = ['\n'];

for r = 1:rows
    
    % add indent for wach row
    for i=1:indent
        txt = [txt '   '];
    end
    
    % add content in each column
    for c = 1:cols
        txt = [txt '\t' msg{r,c}];
    end
    
    % new line for next column
    txt = [txt '\n'];
end

% remove last newline
txt = txt(1:end-2);
end


