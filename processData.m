function dataNI = processData(varargin) %#ok<STOUT>
% processData Process image files from the Penn Natural Image Database,
% calculating texture statistics.
%
%   Options:
%    'forcestats': logical
%       If true, force the regeneration of the statistics file. Otherwise,
%       if the file exists, use it.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('forcestats', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results; %#ok<*NASGU>

% the following variables are used by the scripts below

% names of files to be used are read from imgNamesFile
imgNamesFile='./Natural_Images_Large_NoSky_Index.txt';
% files are located in imgDirectory

imgDirectory='./NaturalImages/';

% similar for the files that will be used to get the whitening filter
imgNamesFile_filter='./Natural_Images_Large_NoSky_Index.txt';
%imgNamesFile_filter='./Natural_Images_TestFilter_Index.txt';
imgDirectory_filter = './NaturalImages/';

% name and directory of output statistics file
analysesFileName='Segmented_Images';
analysesDirectory='./ImageAnalyses/';

% index of image to be considered most in-focus
focusImg = 10;

% pseudocount fraction (?)
PCfrac = 1;

% patch sizes for which to run the analysis
patchSizes = [32, 48];
%patchSizes = [32];
% block averaging factors for which to run the analysis
blockAFs = 2;
%blockAFs = [2, 4];

% collect statistics
% if the statistics files exist already, they are used (unless the
% 'forcestats' option was used)
generateStatFiles

%construct data structure from all analyses & run 2-Gaussian decomposition
loadImageAnalyses
runFocusScript
selectFocusGaussian

% the result structure, dataNI, is defined and updated by the scripts above

end