function safePrint(fname, varargin)
% safePrint Print image to file, avoiding overwrite.
%   safePrint(fname) saves the current figure to a PDF.
%
%   safePrint(fname, 'type', 'png') saves the current figure to a PNG file.
%
%   The extension is added automatically if not provided (which is
%   recommended).
%
%   Options:
%    'printOpts'
%       Options to be passed to Matlab's PRINT.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('type', 'pdf', @(s) ismember(s, {'pdf', 'png'}));
parser.addParameter('printOpts', {}, @(c) iscell(c) && isvector(c));

% show defaults
if nargin == 1 && strcmp(fname, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% add extension if none was given
[~, ~, ext] = fileparts(fname);
if isempty(ext)
    fname = [fname '.' params.type];
end

% if file exists, warn and refuse to overwrite
if exist(fname, 'file')
    warning([mfilename ':fexist'], 'File already exists.');
    return;
end

% set default parameters based on file type
switch params.type
    case 'pdf'
        print('-dpdf', params.printOpts{:}, fname);
    case 'png'
        print('-dpng', '-r300', params.printOpts{:}, fname);
    otherwise
        error([mfilename ':badtype'], 'Unrecognized file type.');
end

end