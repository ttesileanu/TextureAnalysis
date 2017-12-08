function safe_print(fname, type)
% SAFE_PRINT Print image to file, avoiding overwrite.
%   SAFE_PRINT(fname) saves the current figure to a PDF.
%
%   SAFE_PRINT(fname, 'png') saves the current figure to a PNG file.
%
%   The extension is added automatically if not provided (which is
%   recommended).

if nargin < 2
    type = 'pdf';
end

[~, ~, ext] = fileparts(fname);
if isempty(ext)
    fname = [fname '.' type];
end

if exist(fname, 'file')
    warning([mfilename ':fexist'], 'File already exists.');
    return;
end

switch type
    case 'pdf'
        print('-dpdf', fname);
    case 'png'
        print('-dpng', '-r300', fname);
    otherwise
        error([mfilename ':badtype'], 'Unrecognized file type.');
end

end