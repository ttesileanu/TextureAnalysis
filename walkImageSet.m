function res = walkImageSet(fct, imageNames, path, varargin)
% walkImageSet Apply a function to every image in a set and collect the
% results.
%   res = walkImageSet(fct, imageNames, path) loads every image identified
%   by the `imageNames` and located in the folder `path`, and passes it to
%   the function `fct` (after preprocessing by taking the logarithm of each
%   entry; see the options below), together with its index:
%       fct(i, image)
%   The function should return either an empty matrix (which is ignored),
%   or a structure. All the fields in the structures returned for each
%   image are combined into the output `res`.
%
%   The walkImageSet function can take a number of parameters mirroring the
%   parameters of the `preprocessImage` function, which it calls:
%
%   walkImageSet(fct, imageNames, path, masks1, ..., masksN) preprocesses
%   the given masks (which should be cell arrays of matrices) together with
%   the images, and passes them to `fct` as
%       fct(i, image, mask1, ..., maskN)
%
%   walkImgeSet(..., blockAF) downsamples the images by the factor `blockAF
%   before passing them to `fct`. See options below to tweak the
%   downsampling procedure.
%
%   walkImageSet(..., blockAF, filter) filters the images by convolving
%   with fft2(filter) after block averaging. See options below to tweak the
%   filtering procedure.
%
%   walkImageSet(..., blockAF, filter, nLevels) quantizes the images into
%   `nLevels` levels after block-averaging and filtering. This is done by
%   comparing the value of each pixel with the median over the whole image.
%
%   walkImageSet(..., blockAF, filter, nLevels, quantPatchSize) quantizes
%   the image by looking at patches of the given size. See options below to
%   tweak the quantization procedure.
%
%   Options:
%    'averageType': char
%       This is passed to `blockAverage` to set the type of averaging
%       used when downsampling images.
%    'doLog': logical
%       When true, the images are converted to a logarithmic space (this is
%       the default). See `preprocessImage`.
%    'filterType': char
%       Passed to `preprocessImage` to set the type of filtering that is
%       performed (see `filterImage` for a description).
%    'quantType': char
%       Passed to `preprocessImage` to set the type of quantization that is
%       performed (see `quantize` for a description).
%    'threshold': number
%       Threshold to use to avoid taking the logarithm of negative numbers
%       when `doLog` is true. See `preprocessImage` and `convertToLog`.
%    'imageCopies': logical
%       If true, copies of the original, log-transformed, block-averaged,
%       filtered, and final, quantized, images are included in the output
%       structure. The (downsampled) masks are also included. This is true
%       by default for at most 10 images, false otherwise.
%    'progressEvery': double
%       How often to display progress information (in seconds), after the
%       'progressStart' period (see below) elapsed.
%    'progressStart': double
%       How long to wait before displaying progress information for the
%       first time. Set to infinity to never display progress.
%
%   See also: preprocessImage.

% handle the optional masks
firstNonMaskIdx = find(cellfun(@(c) ~iscell(c), varargin), 1);
if isempty(firstNonMaskIdx) || firstNonMaskIdx == 1
    masks = {};
else
    masks = varargin(1:firstNonMaskIdx-1);
    varargin = varargin(firstNonMaskIdx:end);
end

if ~all(cellfun(@length, masks) == length(imageNames))
    error([mfilename ':badmasks'], 'All masks should be cell arrays of the same length as the imageNames.');
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
parser.addOptional('filter', [], @(m) isempty(m) || (ismatrix(m) && isreal(m) && isnumeric(m)));
parser.addOptional('nLevels', checkNumber);
parser.addOptional('quantPatchSize', checkPatchSize);

parser.addParameter('averageType', [], checkStr);
parser.addParameter('doLog', [], checkBool);
parser.addParameter('filterType', [], checkStr);
parser.addParameter('quantType', [], checkStr);
parser.addParameter('threshold', [], checkNumber);

parser.addParameter('imageCopies', [], checkBool);

parser.addParameter('progressEvery', 10, checkNumber);
parser.addParameter('progressStart', 20, checkNumber);

% parse
parser.parse(varargin{:});
params = parser.Results;

% set the default for image copies
if isempty(params.imageCopies)
    params.imageCopies = numel(imageNames) <= 10;
end

t0 = tic;
tEvery = tic;
progressWritten = false;

res = [];
if params.imageCopies
    res.imageCopies.original = cell(size(imageNames));
    res.imageCopies.log = cell(size(imageNames));
    res.imageCopies.blockAveraged = cell(size(imageNames));
    res.imageCopies.filtered = cell(size(imageNames));
    res.imageCopies.final = cell(size(imageNames));
    res.imageCopies.masks = cell(size(imageNames));
end
for i = 1:numel(imageNames)
    % output progress information if required
    if (~progressWritten && toc(t0) > params.progressStart) || ...
            (progressWritten && toc(tEvery) > params.progressEvery)
        disp(['Processing image ' int2str(i) ' of ' int2str(numel(imageNames)) ...
            ', elapsed ' num2str(toc(t0), '%.1f') ' seconds...']);
        progressWritten = true;
        tEvery = tic;
    end
    
    % get the current masks
    crtMasks = cellfun(@(m) m{i}, masks, 'uniform', false);
    
    % preprocess image (and mask, if available)
    crtProcessedMasks = cell(size(crtMasks));
    [crtImage, crtProcessedMasks{:}, ...
        crtOrigImage, crtLogImage, crtAveragedImage, crtFilteredImage] = ...
        preprocessImage(fullfile(path, imageNames{i}), crtMasks{:}, varargin{:});
    
    % call the walker function
    crtRes = fct(i, crtImage, crtProcessedMasks{:});
    
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
    
    % keep track of the source images if asked to do so
    if params.imageCopies
        res.imageCopies.original{i} = crtOrigImage;
        res.imageCopies.log{i} = crtLogImage;
        res.imageCopies.blockAveraged{i} = crtAveragedImage;
        res.imageCopies.filtered{i} = crtFilteredImage;
        res.imageCopies.final{i} = crtImage;
        res.imageCopies.masks{i} = crtMask;
    end
end

if progressWritten
    disp(['Processing finished, took ' num2str(toc(t0), '%.2f') ' seconds.']);
end

res.imageNames = imageNames;
res.path = path;
res.options = copyFields(params, {'averageType', 'blockAF', 'doLog', ...
    'threshold', 'filter', 'filterType', 'quantType', 'quantPatchSize', ...
    'nLevels'});

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