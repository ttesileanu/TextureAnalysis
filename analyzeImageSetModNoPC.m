function [ covM, ev, others ] = analyzeImageSetModNoPC(...
    imgNamesFile, imgDirectory, ...
    imgNamesFile_filter, imgDirectory_filter, ...
    blockAvgFactor, patchSize, ...
    varargin)
% analyzeImageSetModNoPC Split images into patches and calculate texture
% statistics per patch.
%   [covM, ev] = analyzeImageSetModNoPC(
%           imgNamesFile, imgDirectory, imgNamesFile_filter, imgDirectory_filter,
%           blockAvgFactor, patchSize)
%   uses the images identified by 'imgNamesFile_filter' in folder
%   'imgDirectory_filter' to generate a whitening filter (unless the 'filter'
%   option is used -- see below), then uses this filter to whiten and
%   analyze the images identified by 'imgNamesFile' in folder 'imgDirectory'.
%   The analysis splits the images into patches of  size 'patchSize' after
%   having downsampled them by 'blockAvgFactor', and calculates texture
%   statistics for each patch. Each row in the matrix 'ev' contains the
%   values for the 10 independent texture parameters (see getStatistics and
%   processBlock) for each image patch. 'covM' is the covariance matrix for
%   these statistics across the image patches.
%
%   [covM, ev, others] = analyzeImageSetModNoPC(...) also returns
%   information regarding the image position for each patch
%   ('others.imageCoordinates'), the whitening filter ('others.Ffilter'),
%   and per-patch sharpness information ('others.sharpness'). It can also
%   contain the image data that was analyzed to get the texture statistics,
%   provided the 'images' option is set to true. In this case, there will
%   also be fields 'others.origImages', 'others.blockAFImages',
%   'others.whitenedImages', and 'others.binarizedImages'. If the 'images'
%   option is left to its default (which is an empty matrix), the image
%   data will only be kept provided there are at most 25 images to be
%   analyzed.
%
%   analyzeImageSetModNoPC(..., 'segmentation', segmentation) uses the
%   segmentation data from the 'segmentation' structure to focus on
%   particular areas of the images. 'segmentation' needs to have fields
%   'segs', 'segSel', 'segField', and 'segFct'; for a description of the
%   meaning of these fields, see processData.m. When a segmentation is
%   used, the function calculates three sets of statistics: for the full
%   image, for the patches that fall within the segmentation masks, and for
%   the patches that fall outside the segmentation masks.
%
%   [..., details] = analyzeImageSetModNoPc(...) returns details about the
%   images that are being analyzed, like in getImgStats, depending on the
%   options described in the getImgStats help.
%
%   Options:
%    'filter': array
%       Use this matrix as a whitening filter instead of calculating it
%       from images. If this is provided, the values of imgNamesFile_filter
%       and imgDirectory_filter are ignored.
%    'filterFull': bool
%       If true, the whitening filter is applied to the full image. If
%       false, it is applied to each patch separately.
%    'fullImageEv': bool
%       If true, calculate a single set of 10 'ev' values for the whole
%       image (and another two sets 'evF' and 'evB' if a segmentation is
%       present). For now this only works when 'filterFull' is true.
%    'images': bool
%       If true, the function returns 'others.origImages',
%       'others.blockAFImages', 'others.whitenedImages', and
%       'others.binarizedImages' -- cell arrays of the images used in the
%       analyses. If the 'images' option is left to its default (an empty
%       matrix), the image data will only be kept provided there are at
%       most 25 images to be analyzed.
%    'overlappingPatches': bool
%       Set to true to use a pixel-by-pixel sliding patch to evaluate
%       statistics. This generates many more patches but these are no
%       longer independent of each other.
%    'segmentation': struct
%       If this is provided, then the analysis will not only be run on the
%       full images in the dataset, but also on those subsets of the images
%       that are inside the masks provided by the segmentation (see
%       processData.m), and also for the complementary masks. The values
%       will be returned as 'others.evF', 'others.covF', and 'others.evB',
%       'others.covB', respectively.
%    'subSelect': array
%       A list of integers or a boolean mask selecting which images from
%       the analysis (as opposed to filter) set should be processed.
%
%   See also: processBlock, getStatistics, processData, getImgStats.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('segmentation', [], @(s) isempty(s) || isstruct(s));
parser.addParameter('images', [], @(b) isempty(b) || (islogical(b) && isscalar(b)));
parser.addParameter('overlappingPatches', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('subSelect', [], @(v) isempty(v) || (isvector(v) && isnumeric(v)));
parser.addParameter('filter', [], @(m) isempty(m) || (ismatrix(m) && isnumeric(m)));
parser.addParameter('filterFull', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('fullImageEv', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if isempty(params.filter)
    images_filter = parseImageNameFile(imgNamesFile_filter, imgDirectory_filter);
end
imageNames = parseImageNameFile(imgNamesFile, imgDirectory);

if ~isempty(params.subSelect)
    imageNames = imageNames(params.subSelect);
    params.segmentation.segs = params.segmentation.segs(params.subSelect);
    params.segmentation.segSel = params.segmentation.segSel(params.subSelect);
end

if isempty(params.images)
    params.images = numel(imageNames) <= 25;
end

tTotal = tic;
tic

%generate filter (1./sqrt(Fourier Power))
if isempty(params.filter)
    disp('Whitening...');
    [Ffilter] = generateFourierWhitenFilter(imgDirectory_filter, images_filter, ...
        blockAvgFactor, patchSize);
else
    Ffilter = params.filter;
    disp('Using existing filter.');
end

disp(['Whitening took ' num2str(toc, '%.2f') ' seconds.']);

tic;

%compute statistics using only segmented images (a subset of those used to build the filter)
disp('Getting Statistics...');
opts = {...
    'segmentation', params.segmentation, ...
    'images', params.images, ...
    'filterFull', params.filterFull, ...
    'overlappingPatches', params.overlappingPatches, ...
    'fullImageEv', params.fullImageEv};
[ev,numWide,numHigh, others]=getImgStats(imgDirectory, imageNames, blockAvgFactor, ...
    patchSize, Ffilter, opts{:});

if size(ev, 1) > 1
    covM=cov(ev);
else
    % can't calculate (co)variance from single sample
    covM = nan(10);
end
if isfield(others, 'evF')
    if size(ev, 1) > 1
        others.covF = cov(others.evF);
    else
        others.covF = nan(10);
    end
end
if isfield(others, 'evB')
    if size(ev, 1) > 1
        others.covB = cov(others.evB);
    else
        others.covB = nan(10);
    end
end

others.patchStats = struct('numWide', numWide, 'numHigh', numHigh);
others.Ffilter = Ffilter;

disp(['Getting statistics took ' num2str(toc, '%.2f') ' seconds.']);

t = toc(tTotal);
disp([mfilename ' took ' num2str(t, '%.2f') ' seconds to run.']);

end