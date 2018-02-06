function res = walkImageSet(fct, imageNames, path, varargin)
% walkImageSet Apply a function to every image in a set and collect the
% results.
%   res = walkImageSet(fct, imageNames, path) loads every image identified
%   by the `imageNames` and located in the folder `path`, and passes it to
%   the function `fct` (after preprocessing by taking the logarithm of each
%   entry; see the options below), together with its index:
%       fct(i, image, details)
%   The function should return either an empty matrix (which is ignored),
%   or a structure. All the fields in the structures returned for each
%   image are combined into the output `res`. See a description of the
%   `details` structure below.
%
%   res = walkImageSet(fct, imageCount, imageGenerator, ...) uses a function
%   handle (`imageGenerator`) to generate the `imageCount` images that are
%   to be processed. `imageGenerator` is called as `imageGenerator(i)`,
%   where `i` is the index of the image that's being requested, and should
%   return an image that can be handled by `preprocessImage`.
%
%   Furthermore, the walkImageSet function can take a number of parameters
%   mirroring the parameters of the `preprocessImage` function, which it
%   calls:
%
%   walkImageSet(fct, imageNames, path, masks1, ..., masksN) preprocesses
%   the given masks (which should be cell arrays of matrices) together with
%   the images, and passes them to `fct` as
%       fct(i, image, mask1, ..., maskN, details)
%   where `details` is a structure containing `origImage`, `logImage`,
%   `averagedImage`, `filteredImage`, `equalizedImage`, `nonlinearImage`,
%   and `crops`, as returned by `preprocessImage`.
%
%   walkImageSet(..., blockAF) downsamples the images by the factor `blockAF
%   before passing them to `fct`. See options below to tweak the
%   downsampling procedure and other preprocessing options. (The defaults
%   are those from `preprocessImage`.)
%
%   Options:
%    'averageType': char
%       This is passed to `blockAverage` to set the type of averaging
%       used when downsampling images.
%    'doLog': logical
%       When true, the images are converted to a logarithmic space (this is
%       the default). See `preprocessImage`.
%    'equalize': logical or string
%       If set to 'equalize', the image is histogram-equalized after the
%       filtering (and before a potential quantization). If set to
%       'contrast', the image is run through a contrast adaptation
%       algorithm instead (see contrastAdapt.m). Neither is performed if
%       this option is set to `false`.
%    'equalizeType': char
%       This is passed to `equalize` or `contrastAdapt` to set the type of
%       histogram equalization or contrast adaptation that is used.
%    'filter': [], or matrix
%       Whitening filter to use after log and block-averaging, but before
%       equalization and/or quantization. If empty, no whitening is
%       performed.
%    'filterType': char
%       Passed to `preprocessImage` to set the type of filtering that is
%       performed (see `filterImage` for a description).
%    'nonlinearity': [], or vector
%       If non-empty, this instructs preprocessImage to perform a per-pixel
%       nonlinear filtering of the image. This happens after equalization
%       but before quantization. See `applyNonlinearity`.
%    'patchSize': [], int, or [int, int]
%       Patch size to use for `quantize` and/or `equalize`. See the
%       documentation for those two functions for details. If empty, the
%       whole image is processed at once (but see 'quantType' below).
%    'quantize': int
%       If non-empty, perform color quantization to the given number of
%       levels after equalization (if any). The quantization algorithm
%       assumes values in the range [0, 1]; everything outside that range
%       is clipped.
%    'quantizeType': char
%       This can be 'deterministic' (in which case `quantize` is used), or
%       'stochastic' (in which case `stochasticBinarize` is used). Note
%       that 'stochastic' can only be used with binary images, so
%       'quantize' must be 2 in this case.
%    'threshold': number
%       Threshold to use to avoid taking the logarithm of negative numbers
%       when `doLog` is true. See `preprocessImage` and `convertToLog`.
%    'progressEvery': double
%       How often to display progress information (in seconds), after the
%       'progressStart' period (see below) elapsed.
%    'progressStart': double
%       How long to wait before displaying progress information for the
%       first time. Set to infinity to never display progress.
%
%   See also: preprocessImage.

if ~iscell(imageNames)
    if isnumeric(imageNames) && isscalar(imageNames) && isreal(imageNames)
        imageCount = imageNames;
        imageGenerator = path;
        
        imageNames = [];
        path = [];
    else
        error([mfilename ':badimgs'], ['The first argument should either be imageNames, '...
            'a cell array of strings, or imageCount, an integer.']);
    end
else
    imageCount = numel(imageNames);
end

% handle the optional masks
firstNonMaskIdx = find(cellfun(@(c) ~iscell(c), varargin), 1);
if isempty(firstNonMaskIdx) || firstNonMaskIdx == 1
    masks = {};
else
    masks = varargin(1:firstNonMaskIdx-1);
    varargin = varargin(firstNonMaskIdx:end);
end

if ~all(cellfun(@length, masks) == imageCount)
    error([mfilename ':badmasks'], 'All masks should be cell arrays of the same length as the number of images.');
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

checkStr = @(s) isempty(s) || (ischar(s) && isvector(s));
checkBool = @(b) isempty(b) || (islogical(b) && isscalar(b));
checkPatchSize = @(v) isempty(v) || (isnumeric(v) && isvector(v) && ...
    (numel(v) == 1 || numel(v) == 2) && all(v >= 1));
checkNumber = @(x) isempty(x) || (isscalar(x) && isreal(x) && isnumeric(x));

parser.addOptional('blockAF', 1, checkNumber);

parser.addParameter('filter', [], @(m) isempty(m) || (ismatrix(m) && isreal(m) && isnumeric(m)));
parser.addParameter('equalize', 'equalize', @(b) isequal(b, 'false') || ismember(b, {'equalize', 'contrast'}));
parser.addParameter('equalizeType', [], checkStr);
parser.addParameter('patchSize', [], checkPatchSize);
parser.addParameter('averageType', [], checkStr);
parser.addParameter('doLog', [], checkBool);
parser.addParameter('filterType', [], checkStr);
parser.addParameter('nonlinearity', [], @(v) isempty(v) || (isvector(v) && isnumeric(v) && isreal(v)));
parser.addParameter('quantize', [], checkNumber);
parser.addParameter('quantizeType', [], checkStr);
parser.addParameter('threshold', [], checkNumber);

parser.addParameter('progressEvery', 10, checkNumber);
parser.addParameter('progressStart', 20, checkNumber);

% parse
parser.parse(varargin{:});
params = parser.Results;

t0 = tic;
tEvery = tic;
progressWritten = false;

res = [];
preprocessArgs = [{params.blockAF} ...
    structToCell(params, {'averageType', 'doLog', 'equalize', 'equalizeType', ...
    'filter', 'filterType', 'patchSize', 'nonlinearity', 'quantize', 'quantizeType', ...
    'threshold'})];
for i = 1:imageCount
    % output progress information if required
    if (~progressWritten && toc(t0) > params.progressStart) || ...
            (progressWritten && toc(tEvery) > params.progressEvery)
        disp(['Processing image ' int2str(i) ' of ' int2str(imageCount) ...
            ', elapsed ' num2str(toc(t0), '%.1f') ' seconds...']);
        progressWritten = true;
        tEvery = tic;
    end
    
    % get the current masks
    crtMasks = cellfun(@(m) m{i}, masks, 'uniform', false);
    
    % preprocess image (and mask, if available)
    crtProcessedMasks = cell(size(crtMasks));
    if ~isempty(imageNames)
        crtOriginalImage = fullfile(path, imageNames{i});
    else
        crtOriginalImage = imageGenerator(i);
    end
    [crtImage, crtProcessedMasks{:}, crtImages, crtCrops] = ...
        preprocessImage(crtOriginalImage, crtMasks{:}, preprocessArgs{:});
    
    % call the walker function
    details = crtImages;
    details.crops = crtCrops;
    crtRes = fct(i, crtImage, crtProcessedMasks{:}, details);
    
    if ~isempty(crtRes)        
        % add the current stats to the overall structure
        crtFields = fieldnames(crtRes);
        for j = 1:numel(crtFields)
            crtField = crtFields{j};
            crtValue = crtRes.(crtField);
            if isfield(res, crtField)
                res.(crtField) = [res.(crtField) ; crtValue];
            else
                if ischar(crtValue)
                    res.(crtField) = {crtValue};
                else
                    res.(crtField) = crtValue;
                end
            end
        end
    end
end

if progressWritten
    disp(['Processing finished, took ' num2str(toc(t0), '%.2f') ' seconds.']);
end

res.imageNames = imageNames;
res.imageCount = imageCount;
res.path = path;

res.options = copyFields(params, {'averageType', 'blockAF', 'doLog', ...
    'equalize', 'equalizeType', 'filter', 'threshold', 'filterType', ...
    'nonlinearity', 'quantize', 'quantizeType', 'patchSize'});

end

function new_s = copyFields(s, fields)
% copyFields(s, fields) makes a new structure that contains only the chosen
% fields of the original structure. Non-existent fields are silently
% skipped.

new_s = struct;
fields = intersect(fields, fieldnames(s));
for i = 1:length(fields)
    new_s.(fields{i}) = s.(fields{i});
end

end