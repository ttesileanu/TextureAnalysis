function generateLUMFiles(src, dest, varargin)
% generateLUMFiles Convert image files to LUM Matlab files.
%   generateLUMFiles(src, dest) loads all the image files in the source
%   folder 'src' one by one, converts them to a luminance-only format, and
%   stores them into Matlab files with names ending in '_LUM.mat'. Each
%   such file contains a single variable, 'LUM_Image'.
%
%   Currently the function only works with JPEG files whose names end in
%   '.jpg' or '.jpeg' (case-independent).
%
%   Options:
%    'selection' <s/c>
%       Either a cell array of file names to consider in the 'src'
%       directory (as opposed to using 'dir'), or the name of a text file
%       containing the files to focus on. Each line of the text file should
%       have the format 'name=file'; other lines are ignored. Also, lines
%       starting with a '#' are considered comments and ignored.

% parse the optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('selection', [], @(s) (ischar(s) && isvector(s)) || iscell(s));

parser.parse(varargin{:});
params = parser.Results;

% handle reading from file
if ~isempty(params.selection) && ischar(params.selection)
    disp(['Reading file names from ' params.selection '.']);
    f = fopen(params.selection);
    % split into lines
    data = textscan(f, '%s', 'delimiter', '\n');
    fclose(f);
    
    data = data{1}; % focus on the single column
    
    % get rid of insignificant whitespace
    data = cellfun(@strtrim, data, 'uniform', false);
    
    % get rid of comment lines
    data = data(cellfun(@(s) ~isempty(s) && s(1) ~= '#', data));
    
    % find lines with equal signs, and ignore what's to the left
    data = data(cellfun(@(s) sum(s == '=') > 0, data));
    params.selection = cellfun(@(s) s(find(s == '=', 1)+1:end), data, 'uniform', false);
end

if isempty(params.selection)
    % simply use all files in folder
    files = dir(src);
elseif iscell(params.selection)
    files = struct('name', params.selection, 'isdir', false);
end

% some settings
extensions = {'.jpg', '.jpeg'};
postfix = '_LUM.mat';

% loop through the files
for i = 1:length(files)
    crt = files(i);
    
    % skip directories
    if crt.isdir
        disp(['Skipping directory ' crt.name '...']);
        continue;
    end
    
    [~, fname, ext] = fileparts(crt.name);
    if any(strcmpi(ext, extensions))
        % this is an image file
        LUM_Image = imageToLUM(fullfile(src, crt.name));
        disp(['Converted ' crt.name ', size ' int2str(size(LUM_Image, 2)) ' x ' ...
            int2str(size(LUM_Image, 1)) '...']);
        save(fullfile(dest, [fname postfix]), 'LUM_Image');
    else
        disp(['Skipping unrecognized file ' crt.name '...']);
    end
end

end