function images = parseImageNameFile(imgNameFile, path)
% parseImageNameFile Parse file of image names.
%   images = parseImageNameFile(imgNameFile) parses a text file containing
%   image names (one per line) and returns a cell array containing all the
%   names. Lines starting with the pound sign '#' are ignored.
%
%   parseImageNameFile(imgNameFile, path) concatenates the given path in
%   front of each name found in the file. This uses Matlab's FULLFILE
%   command, so it automatically adds the appropriate file separator
%   character when necessary.
%
%   See also: FULLFILE.

% handle both forms of the command
if nargin < 2
    path = '';
end

% read the file
f = fopen(imgNameFile);
raw = textscan(f, '%s', 'delimiter', '\n\r');
fclose(f);

% filter out comment lines
raw = raw{1};

n = 1;
images = cell(1, length(raw));
for i = 1:length(raw)
    if raw{i}(1) == '#'
        % comment line
        continue;
    end
    images{n} = fullfile(path, strtrim(raw{i}));
    n = n + 1;
end

% we might have preallocated too much if there are comment lines
images = images(1:n-1);

end