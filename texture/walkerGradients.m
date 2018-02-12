function crtRes = walkerGradients(i, image, varargin)
% walkerGradients A walker function for walkImageSet that calculates
% summary statistics of gradient filters.
%   crtRes = walkerGradients(i, image, [mask1, ..., maskN], details) runs
%   `analyzePatchGradients` on the given preprocessed image and outputs the
%   relevant statistics. This function is meant for use with `walkImageSet`.
%
%   See also: `walkImageSet`, `analyzePatchGradients`.

% check for images preprocessed out of existence
if isempty(image)
    % if preprocessing reduced the image to nothing, skip it
    warning([mfilename ':smallimg'], ['Image #' int2str(i) ' was trimmed to zero size by preprocessing -- skipping.']);
    crtRes = [];
    return;
end

% pick out masks and details
details = varargin{end};
masks = varargin(1:end-1);

if isempty(masks)
    % need to create a mask as large as the unscaled & untrimmed
    % image
    mask = ones(details.crops.final(3), details.crops.final(4));
else
    if length(masks) > 1
        error([mfilename ':nmasks'], 'Can use at most one mask.');
    end
    mask = masks{1};
    % if the current mask is empty, skip image
    if isempty(mask)
        crtRes = [];
        return;
    end
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('analysisPatchSize', [], @(x) isnumeric(x) && (isscalar(x) || ...
    (isvector(x) && length(x) == 2)));

% parse
parser.parse(details.args{:});
params = parser.Results;

% calculate statistics
try
    crtRes = analyzePatchGradients(image, params.analysisPatchSize, mask, ...
        'maskCrop', details.crops.final);
catch errMsg
    newErrMsg.message = ['Image #' int2str(i) ': ' errMsg.message];
    newErrMsg.identifier = errMsg.identifier;
    newErrMsg.stack = errMsg.stack;
    error(newErrMsg);
end
if isempty(crtRes.ev)
    crtRes = [];
    return;
end
crtRes.imgIds = i*ones(size(crtRes.ev, 1), 1);

crtRes = rmfield(crtRes, intersect(...
    {'patchSize', 'overlapping', 'minPatchUsed', 'maxPatchesPerImage'}, ...
    fieldnames(crtRes)));
if i > 1
    crtRes = rmfield(crtRes, intersect({'statsDesc'}, fieldnames(crtRes)));
end

% keep track of the source images if asked to do so
% if details.params.imageCopies
%     crtRes.imageCopies.original = details.origImage;
%     crtRes.imageCopies.log = details.logImage;
%     crtRes.imageCopies.blockAveraged = details.averagedImage;
%     crtRes.imageCopies.filtered = details.filteredImage;
%     crtRes.imageCopies.equalized = details.equalizedImage;
%     crtRes.imageCopies.nonlinear = details.nonlinearImage;
%     crtRes.imageCopies.final = image;
%     crtRes.imageCopies.mask = mask;
%     crtRes.imageCopies.crops = details.crops;
% end

end