function [power, radii, fourier] = powerspec(im, varargin)
% POWERSPEC Calculate image power spectrum.
%   [power, radii] = POWERSPEC(im) calculates the radial power spectrum for
%   the image. If the input image is a 3d array, the power spectrum is
%   calculated for each slice through the third dimension, then averaged.
%
%   POWERSPEC(im, 'nbins', nbins) chooses the number of bins to use for the power
%   spectrum.
%
%   [..., fourier] = POWERSPEC(im) also returns the 2d Fourier transform,
%   shifted so that the zero is in the center.
%
%   POWERSPEC(fourier, 'fromFourier', true) calculates the radial power
%   spectrum given the 2d Fourier transform.
%
%   Options:
%    'nbins'
%       Number of bins to use to calculate radial power spectrum.
%    'patchSize'
%       Patch size, as a single scalar or a pair [patchSizeY, patchSizeX]
%       (note that y size is first!). If given, the power spectrum is
%       calculated for each patch and averaged.
%    'fromFourier'
%       If set to `true`, the first argument is assumed to be a 2d Fourier
%       transform of the image. In this case, `patchSize` is ignored.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('nbins', 100, @(n) isnumeric(n) && isscalar(n) && n > 1);
parser.addParameter('fromFourier', false, @(b) islogical(b) && isscalar(b));

% show defaults
if nargin == 1 && strcmp(im, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

nPatches = size(im, 3);

% get the average Fourier spectrum (or use the provided one)
if ~params.fromFourier
    fourier = zeros(size(im, 1), size(im, 2));
    for i = 1:nPatches
        fourier = fourier + fftshift(fft2(im(:, :, i)));
    end
    fourier = fourier / nPatches;
else
    fourier = im;
end

% calculate Fourier power
fourierAbs2 = abs(fourier).^2;

% find distances from center to each of the pixels in the Fourier power matrix
xMax = size(im, 2)/2;
yMax = size(im, 1)/2;
xgv = -xMax+1:xMax;
ygv = -yMax+1:yMax;
[x, y] = meshgrid(xgv, ygv);

r2 = x.^2 + y.^2;
r1 = sqrt(r2);
maxRadius = min(xMax, yMax);
radii = linspace(0, maxRadius, params.nbins);

% collect the power in each ring
power0 = accumarray(1 + floor(r1(:)*length(radii)/maxRadius), fourierAbs2(:));
counts = accumarray(1 + floor(r1(:)*length(radii)/maxRadius), ones(numel(fourierAbs2), 1));
power = power0 ./ counts;
power = power(1:length(radii));

end