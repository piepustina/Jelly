function anonymousFunction(varargin)
% This will write an anonymous function a an m file.
% Note that this will only work if the function does not have capture
% variables. So this will not work:
%   a = 100;
%   func = @(x) x+a;
%   write_anon_function_to_file(func)
% But this will:
%   func = @(x) x+100;
%   write_anon_function_to_file(func)
% Determine the function name.

switch length(varargin)
    case 2
        [filepath, name, ~] = fileparts(varargin{2});
    otherwise
        name = 'myfunc';
end

% Extract and parse the function contents.
txt = formattedDisplayText(varargin{1}); % capture the output
txt = char(txt); % convert to char
txt = strtrim(txt); % remove leading and trailing whitespace
ind = strfind(txt,')'); ind = ind(1); % find the end of the input arguments

function_body = txt((ind+1):end);

% Compose the function and write to file.
txt = sprintf('function out = %s(%s)\nout = %s;\nend', ...
    name,txt(3:(ind-1)), function_body);
writelines(txt, string(filepath) + "/" + string(name) + ".m");


end