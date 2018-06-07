function [power, radii, fourier] = powerspec(im, varargin)
% POWERSPEC Calculate image power spectrum.
%   [power, radii] = POWERSPEC(im) calculates the radial power spectrum for
%   the image.
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
parser.addParameter('patchSize', [], @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) <= 2));
parser.addParameter('fromFourier', false, @(b) islogical(b) && isscalar(b));

if nargin == 1 && strcmp(im, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle scalar patch size and no patch size
if ~isempty(params.patchSize) && isscalar(params.patchSize)
    params.patchSize = [params.patchSize params.patchSize];
elseif isempty(params.patchSize)
    % no patch size is equivalent to one image-sized patch
    params.patchSize = size(im);
end

if ~isempty(params.patchSize)
    npatches = floor(size(im) ./ params.patchSize);
end

if ~params.fromFourier
    fourier = zeros(params.patchSize);
    % fourier_abs2 = zeros(params.patchSize);
    for i = 1:npatches(1)
        for j = 1:npatches(2)
            crt_patch = im(1+(i-1)*params.patchSize(1):i*params.patchSize(1), ...
                1+(j-1)*params.patchSize(2):j*params.patchSize(2));
            crt_fourier = fftshift(fft2(crt_patch));
            %         crt_fourier_abs2 = abs(crt_fourier).^2;
            
            fourier = fourier + crt_fourier;
            %         fourier_abs2 = fourier_abs2 + crt_fourier_abs2;
        end
    end
    fourier = fourier / prod(npatches);
else
    fourier = im;
end

fourier_abs2 = abs(fourier).^2;

xmax = params.patchSize(2)/2;
ymax = params.patchSize(1)/2;
xgv = -xmax+1:xmax;
ygv = -ymax+1:ymax;
[x, y] = meshgrid(xgv, ygv);

r2 = x.^2 + y.^2;
r1 = sqrt(r2);
max_radius = min(xmax, ymax);
radii = linspace(0, max_radius, params.nbins);

power0 = accumarray(1 + floor(r1(:)*length(radii)/max_radius), fourier_abs2(:));
counts = accumarray(1 + floor(r1(:)*length(radii)/max_radius), ones(numel(fourier_abs2), 1));
power = power0 ./ counts;
power = power(1:length(radii));

end