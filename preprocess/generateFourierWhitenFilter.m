function Ffilter = generateFourierWhitenFilter(images, patchSize, varargin)
% generateFourierWhitenFilter Generate a filter to remove average two-point
% correlations from images.
%   Ffilter = generateFourierWhitenFilter(images, patchSize) generates a
%   filter of the given size that can be used to remove the average
%   two-point correlations from images. The `patchSize` can be either a
%   single integer or a pair of integers to use non-square patches (note
%   that the order is [rows, columns]). Before calculating the
%   autocorrelation spectrum, the images are preprocessed according to the
%   'preprocessing' option (see below). Any other key-value arguments are
%   directly passed to walkImages.
%
%   The returned matrix, Ffilter, is the Fourier transform of the filter
%   that must be used in a convolution to whiten the image. It is
%   calculated as the matrix of inverse square roots of the Fourier power
%   values averaged over all the patches in the image set.
%
%   Options:
%    'preprocessing':
%       A cell array of preprocessing functions that will be applied to
%       each image, in the same format as for walkImages.
%    Any other key-value arguments are directly passed to walkImages.
%
%   See also: walkImages.

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

% handle square or rectangular patch sizes
if numel(patchSize) == 1
    patchSize = [patchSize patchSize];
end

% walk through the image set and sum the power spectra
Fpower = zeros(patchSize);
count = 0;

walkImages([params.preprocessing {@walker}], images, unmatched);

% take the average
Fpower = Fpower / count;
% this is the Fourier transform of the filter that removes the average autocorrelation
Ffilter = 1 ./ sqrt(Fpower);

    function result = walker(~, image, ~)
        % process each image, sum power spectra in each patch
        patchifier = ImagePatchifier(image, patchSize);
        
        while patchifier.next
            patch = patchifier.get;
            patchF = fft2(patch);
            Fpower = Fpower + (patchF .* conj(patchF));
        end
        count = count + patchifier.numPatches;
                
        result = [];
    end

end