function safe_print(fname, varargin)
% SAFE_PRINT Print image to file, avoiding overwrite.
%   SAFE_PRINT(fname) saves the current figure to a PDF.
%
%   SAFE_PRINT(fname, 'type', 'png') saves the current figure to a PNG file.
%
%   The extension is added automatically if not provided (which is
%   recommended).
%
%   Options:
%    'printopts'
%       Options to be passed to Matlab's PRINT.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('type', 'pdf', @(s) ismember(s, {'pdf', 'png'}));
parser.addParameter('printopts', {}, @(c) iscell(c) && isvector(c));

if nargin == 1 && strcmp(fname, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

[~, ~, ext] = fileparts(fname);
if isempty(ext)
    fname = [fname '.' params.type];
end

if exist(fname, 'file')
    warning([mfilename ':fexist'], 'File already exists.');
    return;
end

switch params.type
    case 'pdf'
        print('-dpdf', params.printopts{:}, fname);
    case 'png'
        print('-dpng', '-r300', params.printopts{:}, fname);
    otherwise
        error([mfilename ':badtype'], 'Unrecognized file type.');
end

end