function dataNI = processData_Penn(varargin)
% processData_Penn Process image files from the Penn Natural Image Database,
% calculating texture statistics.
%
%   This calls 'processData' with appropriate options; any options passed
%   to this function are redirected to 'processData'.
%
%   See also: processData.

% names of files to be used are read from imgNamesFile
imgNamesFile='./Natural_Images_Large_NoSky_Index.txt';
% files are located in imgDirectory
imgDirectory='NaturalImages/';

% similar for the files that will be used to get the whitening filter
%imgNamesFile_filter='./Natural_Images_Large_NoSky_Index.txt';
%imgDirectory_filter = './NaturalImages/';
imgNamesFile_filter='./Natural_Images_Large_NoSky_Index.txt';
imgDirectory_filter='NaturalImages/';

% name and directory of output statistics file
analysesDirectory='./PennAnalyses/';

% which image to be considered most in-focus
focusImg = 10;

dataNI = processData(imgNamesFile, imgDirectory, ...
    imgNamesFile_filter, imgDirectory_filter, ...
    analysesDirectory, 'focusImg', focusImg, ...
    'patchSizes', [32 48], 'blockAFs', 2, ...
    varargin{:});

end