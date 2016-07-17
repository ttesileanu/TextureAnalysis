function [image, origImage] = preprocessImage(image0, blockAF, varargin)
% preprocessImage Preprocess the image by going to log space and block
% averaging.
%   image = preprocessImage(image0, blockAF) preprocesses the image by
%   taking the logarithm of the entries and block averaging using the given
%   factor. 'image0' can either be the filename for an image, or a matrix
%   containing the image data itself.
%
%   [image, origImage] = preprocessImage(...) also returns the original
%   image, before the logarithm and the block averaging. This is useful if
%   the input image is given as a file name, and the caller wants to have
%   access to the original contents of the file.
%
%   Options:
%    'threshold': double
%       Entries below this number of the input image are replaced by this
%       threshold before taking the logarithm.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('threshold', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);

% parse
parser.parse(varargin{:});
params = parser.Results;

if ischar(image0)
    image = loadLUMImage(image0);
else
    image = image0;
end

origImage = image;

%Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
%LUM_Image(LUM_Image <= 0) = 10^(-.356);
image(image < params.threshold) = params.threshold;
%image(image <= 0) = params.threshold;

image = log(image);
image = blockAverage(image, blockAF, 'avg');

end