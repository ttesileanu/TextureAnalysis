function maskNew = fixContourVoc(mask, varargin)
% fixContourVoc Assign contours in VOC segmentations to nearby objects.
%   maskNew = fixContourVoc(mask) assigns contour pixels to nearby objects
%   in segmentations from the VOC database. For each contour pixel, a patch
%   of size 10 (see options below to change this) is analyzed to find the
%   most salient nearby object. The contour pixels are replaced with this
%   object.
%
%   Contours are identified by having value 255. Value 0 is considered
%   background.
%
%   Options:
%    'patchSize'
%       The size of the patch to use when searching for nearby objects.
%       (default: 10)
%    'backgroundAsObject'
%       Set to false to avoid treating the background (value 0) as just
%       another color.
%       (default: true)

if isempty(mask)
    maskNew = mask;
    return;
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('patchSize', 10, @(x) isscalar(x) && isnumeric(x) && isreal(x) && x > 1);
parser.addParameter('backgroundAsObject', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% find locations of contour pixels
[cRows, cCols] = ind2sub(size(mask), find(mask == 255));

maxObjId = max(mask(mask < 255));
maskNew = mask;
hs = floor(params.patchSize/2);
for i = 1:length(cRows)
    row = cRows(i);
    col = cCols(i);
    
    rows = max(row-hs, 1):min(row+hs, size(mask, 1));
    cols = max(col-hs, 1):min(col+hs, size(mask, 2));
    
    patch = mask(rows, cols);
    counts = histc(patch(:), ((1-params.backgroundAsObject):(maxObjId+1)) - 0.5);
    counts = counts(1:end-1);
    if any(counts > 0)
        [~, newId] = max(counts);
        maskNew(row, col) = newId - params.backgroundAsObject;
    else
        maskNew(row, col) = 0;
    end
end

end