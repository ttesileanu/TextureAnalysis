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

parser.addParamValue('selection', [], @(s) (ischar(s) && isvector(s)) || iscell(s));

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
        image = double(imread(fullfile(src, crt.name)));
        
        disp(['Converting ' crt.name ', size ' int2str(size(image, 2)) ' x ' ...
            int2str(size(image, 1)) '...']);
        if size(image, 3) == 3
            % convert to gray scale
            % this is a simple transformation that approximates the proper
            % gamma correction for sRGB
            image = sqrt(0.299*image(:, :, 1).^2 + ...
                         0.587*image(:, :, 2).^2 + ...
                         0.114*image(:, :, 3).^2);
        elseif size(image, 3) ~= 1
            % don't really know what to do with other numbers of samples
            error([mfilename ':badsamps'], 'Can only work with 1 or 3 samples per pixel.');
        end
        
        % store this to a Matlab file
        LUM_Image = image; %#ok<NASGU>
        save(fullfile(dest, [fname postfix]), 'LUM_Image');
    else
        disp(['Skipping unrecognized file ' crt.name '...']);
    end
end

end