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
%   Use the 'labelsDirectory' option below to read in labels for the
%   segmented objects.
%
%   Options:
%    'labelsDirectory': string
%       Directory where to search for XML files describing the objects in
%       the segmentations.
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
parser.addParameter('labelsDirectory', '', @(s) isempty(s) || (ischar(s) && isvector(s)));

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
        
        % is there any label data?
        if ~isempty(params.labelsDirectory)
            labelsFile = fullfile(params.labelsDirectory, [file_root '.xml']);
            xmlData = [];
            try
                xmlData = parseXML(labelsFile);
            catch
                 warning([mfilename ':noxml'], ['XML data could not be loaded for ' segs(i).image '.']);   
            end
            if ~isempty(xmlData)
                [FG.objMatLabels, FG.objMatDetails] = translateXml(xmlData);
            end
        end
    end
    
    segs(i).FG = FG;
end

end

function [labels, details] = translateXml(xml)
% Translate the XML structure into object labels and other details that may
% be available in the file (bounding box, pose, etc.)

labels = {};
details = {};

if ~strcmp(xml.Name, 'annotation')
    return;
end

for i = 1:numel(xml.Children)
    node = xml.Children(i);
    if strcmp(node.Name, 'object')
        objDetails = struct;
        for j = 1:numel(node.Children)
            subnode = node.Children(j);
            if strcmp(subnode.Name, 'name')
                % assuming this is properly formatted, and keeping only
                % information from the last 'name' tag from every 'object'
                objName = subnode.Children(1).Data;
            end
            if strcmp(subnode.Name, 'pose')
                objDetails.pose = subnode.Children(1).Data;
            end
            if any(strcmp(subnode.Name, {'truncated', 'occluded', 'difficult'}))
                objDetails.(subnode.Name) = logical(subnode.Children(1).Data);
            end
            if strcmp(subnode.Name, 'bndbox')
                bndBox = struct;
                for k = 1:numel(subnode.Children)
                    subsubnode = subnode.Children(k);
                    if ~strcmp(subsubnode.Name, '#text')
                        bndBox.(subsubnode.Name) = str2double(subsubnode.Children(1).Data);
                    end
                end
                objDetails.bndBox = bndBox;
            end
        end
        labels = [labels objName]; %#ok<AGROW>
        details = [details objDetails]; %#ok<AGROW>
    end
end

end