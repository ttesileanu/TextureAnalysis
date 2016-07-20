function dataNI = processData(...
    imgNamesFile, imgDirectory, ...
    imgNamesFile_filter, imgDirectory_filter, analysesDirectory, ...
    varargin)
% processData Process image files, calculating texture statistics, and save
% results to file.
%   dataNI = processData(imgNamesFile, imgDirectory,
%                        imgNamesFile_filter, imgDirectory_filter,
%                        analysesDirectory)
%   returns a struct containing texture analyses for the files listed in
%   'imgNamesFile' and located in 'imgDirectory', using a whitening filter
%   based on the files listed in 'imgNamesFile_filter' and located in
%   'imgDirecory_filter'. The statistics are also stored in a file whose
%   name starts with 'stats', inside the 'analysesDirecotry' (the prefix
%   for this file can be changed; see options below). If the statistics
%   files already exist, an error is generated, unless 'forceSave' is true.
%   The whitening filter from the existing files can be used (if available)
%   if 'filters' is set to 'load'.
%
%   processData(..., 'segs', segs) uses the segmentation given in the
%   structure array 'segs'. 'segs' can also be the name of a .mat file
%   containing the 'segs' (and optionally 'segSel' or 'seg_select') data.
%   Within the 'segs' structure, the data from field 'fgMat' is used, unless
%   a different one is chosen (see options below). If instead the structure
%   'segs' has a field 'FG', then 'fgMat' is searched for within 'FG'. Each
%   entry in segs can also be a structure array (as opposed to a single
%   structure); in this case, the first element in the array is used, unless
%   the appropriate option below is used. The segmentation data can be
%   processed before using by employing the 'segFct' option.
%
%   Options:
%    'blockAFs': vector of int
%       Block averaging sizes used in the analysis.
%       (default: 1)
%    'filters': 'load', or cell array of arrays
%       If equal to 'load' (default), the function attempts to load the
%       whitening filters from old analysis files -- if these exist. This
%       option can also be used to directly provide the filters, in the
%       form of a cell array of matrices. In both of these cases,
%       imgNamesFile_filter and imgDirectory_filter are ignored.
%    'filterFull': bool
%       If true, the whitening filter and the binarization are applied to
%       the full image. If false, they are applied to each patch separately.
%       NOTE: full-image binarization is very slow (over 100x slower than
%       per-patch binarization).
%    'focusImg': int
%       The index of an image whose patches are mostly in focus. This is
%       used with the gaussian mixture fit to find the component
%       corresponding to in-focus patches.
%       (default: 1)
%    'forceSave': logical
%       If true, force saving the analysis results. Otherwise the function
%       refuses to overwrite existing analyses.
%    'fullImageEv': bool
%       If true, calculate a single set of 10 'ev' values for the whole
%       image (and another two sets 'evF' and 'evB' if a segmentation is
%       present). For now this only works when 'filterFull' is true.
%    'images': bool
%       If true, include the grayscale images, both block-averaged and
%       original versions, plus the whitened image, and the binarized image
%       in the output structure. This is true by default if there are at
%       most 25 images, and is false by default otherwise.
%    'loadStats': bool
%       If true, load statistics from file (with the name implied by the
%       'outPrefix' option), provided the file exists and has compatible
%       structure. Note that this can result in wrong results because it's
%       impossible to check that nothing changed between the time when the
%       file was generated and the time when it is used. By default, the
%       processData only uses the whitening filters from file but
%       recalculates the statistics.
%    'outPrefix': char
%       Prefix for output statistics file.
%       (default: 'stats' if 'fullFilter' is false, 'stats_fullflt' else)
%    'overlappingPatches': bool
%       Set to true to use a pixel-by-pixel sliding patch to evaluate
%       statistics. This generates many more patches but these are no
%       longer independent of each other.
%    'patchSizes': vector of int
%       Patch size(s) to be used for the analysis.
%       (default: 32)
%    'segs': string, or segmentation structure
%       Either name of file where segmentation data is stored, or a
%       segmentation structure, as described above.
%    'segFct': function handle
%       A function to be applied to the 'segField' segmentation data before
%       use.
%    'segField': char
%       Name of field to use from the 'segs' structure.
%       (default: 'fgMat')
%    'segSel': array
%       An array identifying which segmentation to use for each image. A
%       value of NaN results in skipping the image.
%       (default: use first segmentation for all images)
%    'subSelect': array
%       List of images to use, or a boolean mask indicating which images to
%       be kept. The determination of the focus gaussian doesn't run when
%       this list is nonempty.
%    'progressEvery': float
%       How often to display progress information (in seconds), after the
%       'progressStart' period (see below) elapsed.
%    'progressStart': float
%       How long to wait before displaying progress information for the
%       first time. Set to infinity to never display progress.
%
%   See also: analyzeImageSetModNoPC, getImgStats.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('blockAFs', 1, @(v) isvector(v) && isnumeric(v));
parser.addParameter('filters', 'load', @(c) isempty(c) || iscell(c) || ...
    strcmp(c, 'load'));
parser.addParameter('filterFull', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('focusImg', 1, @(i) isscalar(i) && isnumeric(i));
parser.addParameter('forceSave', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('fullImageEv', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('images', [], @(b) isempty(b) || (islogical(b) && isscalar(b)));
parser.addParameter('loadStats', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('outPrefix', [], @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParameter('overlappingPatches', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('patchSizes', 32, @(v) isvector(v) && isnumeric(v));
parser.addParameter('segs', [], @(s) isempty(s) || (ischar(s) && isvector(s)) || ...
    (isstruct(s) && isvector(s)));
parser.addParameter('segFct', [], @(f) isempty(f) || isa(f, 'function_handle'));
parser.addParameter('segField', 'fgMat', @(s) ischar(s) && isscalar(s));
parser.addParameter('segSel', [], @(v) isempty(v) || (isvector(v) && isnumeric(v)));
parser.addParameter('subSelect', [], @(v) isempty(v) || (isvector(v) && ...
    (islogical(v) || isnumeric(v))));
parser.addParameter('progressEvery', 10, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('progressStart', 20, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% these will be used by 'generateStatFiles'
params.imgNamesFile = imgNamesFile;
params.imgDirectory = imgDirectory;

params.imgNamesFile_filter = imgNamesFile_filter;
params.imgDirectory_filter = imgDirectory_filter;

params.analysesDirectory = analysesDirectory;

if isempty(params.outPrefix)
    if params.filterFull
        params.outPrefix = 'stats_fullflt';
    else
        params.outPrefix = 'stats';
    end
end

if ischar(params.segs)
    [~, params.segName, ~] = fileparts(params.segs);
    segData = open(params.segs);
    params.segs = segData.segs;
    % use 'segSel' from file only if it wasn't overridden by an option
    if isempty(params.segSel)
        if isfield(segData, 'segSel')
            params.segSel = segData.segSel;
        elseif isfield(segData, 'seg_select')
            params.segSel = segData.seg_select;
        end
    end
else
    params.segName = '';
end

% collect statistics
dataNI = generateStatFiles(params);

% change some field names in dataNI and exchange some entries of the
% covariance matrix to match conventions from psychophysics data
dataNI = convertData(dataNI);

% run 2-Gaussian decomposition to find blurred and in-focus patches
if isempty(params.subSelect)
    dataNI = runFocusScript(dataNI, 'focusImage', params.focusImg);
end

end