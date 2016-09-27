function Ffilter = generateFourierWhitenFilter(imageNames, path, ...
    blockAvgFactor, patchSize, varargin)
% generateFourierWhitenFilter Generate a whitening filter.
%   Ffilter = generateFourierWhitenFilter(imageNames, path, ...
%               blockAvgFactor, patchSize)
%   generates a whitening filter using the images with the given
%   `imageNames` found in the given `path`. The images are first downsampled
%   using `blockAvgFactor`, then split into patches of (linear) size given
%   by `patchSize`. `patchSize` can be either a single integer or a pair of
%   integers to use non-square patches. In this case, the order is row (y)
%   size first, and then column (x) size, `[patchSizeY, patchSizeX]`!
%
%   The returned matrix, Ffilter, is the matrix of inverse square roots of
%   the Fourier power values averaged over all the patches in the image
%   set -- i.e., it is the 2d Fourier transform of the filter that must be
%   used in a convolution to whiten an image.
%
%   Options:
%    'progressEvery': float
%       How often to display progress information (in seconds), after the
%       'progressStart' period (see below) elapsed.
%    'progressStart': float
%       How long to wait before displaying progress information for the
%       first time. Set to infinity to never display progress.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('progressEvery', [], @(x) isnumeric(x) && isscalar(x));
parser.addParameter('progressStart', [], @(x) isnumeric(x) && isscalar(x));

% parse -- just for error checking, since we directly pass the arguments to
% walkImageSet, anyway
parser.parse(varargin{:});

if numel(patchSize) == 1
    patchSize = [patchSize patchSize];
end

Fpower = zeros(patchSize);
count = 0;

walkImageSet(@walker, imageNames, path, blockAvgFactor, varargin{:});

Fpower = Fpower / count;
Ffilter = 1 ./ sqrt(Fpower);

    function res = walker(~, image)
        nPatches = floor(size(image) ./ patchSize);
        
        for k = 1:nPatches(1)
            rows = ((k-1)*patchSize(1) + 1):(k*patchSize(1));
            for l = 1:nPatches(2)
                cols = ((l-1)*patchSize(2) + 1):(l*patchSize(2));
                patchF = fft2(image(rows, cols));
                Fpower = Fpower + (patchF .* conj(patchF));
            end
        end
        count = count + prod(nPatches);
        
        res = [];
    end

end