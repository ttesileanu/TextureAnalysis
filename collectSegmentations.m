function [segs, segSel] = collectSegmentations(imgNamesFile, segDirectory, varargin)
% collectSegmentations Collect image segmentations into a single structure.
%   segs = collectSegmentations(imgNamesFile, segDirectory) loads
%   segmentation data from 'segDirectory' corresponding to the files from
%   the 'imgNamesFile'. Matching segmentation files are assumed to have the
%   same name (but different extension). When a segmentation file is not
%   found, a warning is issued, and a blank segmentation is returned.
%
%   [segs, segSel] = collectSegmentations(...) returns a selection vector
%   for the segmentations in which all entries are equal to 1 except when a
%   segmentation file is not found, in which case the entry is set to NaN.
%
%   Options:
%    'progressEvery': float
%       How often to display progress information (in seconds), after the
%       'progressStart' period (see below) elapsed.
%    'progressStart': float
%       How long to wait before displaying progress information for the
%       first time. Set to infinity to never display progress.
%    'warnNotFound': bool
%       If true, issue a warning for segmentations that aren't found.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('progressEvery', 10, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('progressStart', 20, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('warnNotFound', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

images = parseImageNameFile(imgNamesFile);

segs = repmat(struct, length(images), 1);
segSel = ones(length(images), 1);

t0 = tic;
tLast = tic;

for i = 1:length(images)
    if toc(t0) > params.progressStart
        if toc(tLast) > params.progressEvery
            disp([mfilename ' working on image ' int2str(i) '/' int2str(length(images)) ...
                ', ' images(i).path ...
                '... ' num2str(toc(t0), '%.2f') ' seconds elapsed.']);
            tLast = tic;
        end
    end
    
    segs(i).image = images(i).path;
    
    % search for files with same name
    FG = struct;
    [~, file_root] = fileparts(segs(i).image);
    possible_files = dir(fullfile(segDirectory, [file_root '.*']));
    
    % skip directories
    possible_files = possible_files(~[possible_files.isdir]);
    
    % XXX could also filter by extension, keeping only PNG/GIF/TIFF/...
    
    if isempty(possible_files)
        if params.warnNotFound
            warning([mfilename ':notfound'], ['Segmentation file for ' segs(i).image ' not found.']);
        end
        FG.segPath = '';
        FG.fgMat = [];
        FG.contourMat = [];
        FG.objMat = [];
        segSel(i) = nan;
    else
        FG.segPath = possible_files(1).name;
        if numel(possible_files) > 1
            warning([mfilename ':notfound'], ['Several segmentation files found for ' segs(i).image ...
                ' Using ' FG.segPath '.']);
        end
        segData = imread(fullfile(segDirectory, FG.segPath));
        
        FG.objMat = segData;
        FG.contourMat = (segData == 255);
%        FG.fgMat = (segData > 0 & segData < 255);
        % include contour in the foreground
        FG.fgMat = (segData > 0);
    end
    
    segs(i).FG = FG;
end

end