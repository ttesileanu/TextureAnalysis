function Ffilter = generateFourierWhitenFilter(images, varargin)
% generateFourierWhitenFilter Generate a filter to remove average two-point
% correlations from images.
%   Ffilter = generateFourierWhitenFilter(images) generates a filter that
%   can be used to remove the average two-point correlations from images.
%   Before calculating the autocorrelation spectrum, the images are
%   preprocessed according to the 'preprocessing' option (see below). The
%   filter size is determined by the size of the first two dimensions of
%   the image returned from the last preprocessing function. Use the
%   `patchify` function during preprocessing to split input images into
%   patches of equal sizes.
%
%   The returned matrix, Ffilter, is the Fourier transform of the filter
%   that must be used in a convolution to whiten the image (this is the
%   form that can be used directly with `filterImage`). It is calculated as
%   the matrix of inverse square roots of the Fourier power values averaged
%   over all the patches in the image set.
%
%   Options:
%    'preprocessing':
%       A cell array of preprocessing functions that will be applied to
%       each image, in the same format as for `walkImages`.
%    Any other key-value arguments are directly passed to `walkImages`.
%
%   See also: walkImages, patchify, filterImage.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('preprocessing', {}, @(c) iscell(c) && (isempty(c) || isvector(c)));

% show defaults
if nargin == 1 && strcmp(images, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% walk through the image set and sum the power spectra
fPower = [];
count = 0;

walkImages([params.preprocessing {@walker}], images, unmatched);

% take the average
fPower = fPower / count;
% this is the Fourier transform of the filter that removes the average autocorrelation
Ffilter = 1 ./ sqrt(fPower);

    function result = walker(~, image, ~)
        % process each image patch, sum power spectra in each patch
        nPatches = size(image, 3);
        imPower = zeros(size(image, 1), size(image, 2));
        for i = 1:nPatches
            patchF = fft2(image(:, :, i));
            imPower = imPower + patchF .* conj(patchF);
        end

        % combine with overall sum for all images
        if isempty(fPower)
            fPower = imPower;
        else
            fPower = fPower + imPower;
        end
        count = count + nPatches;
                
        result = [];
    end

end