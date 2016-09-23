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

parser.addParameter('progressEvery', 10, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('progressStart', 20, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results; %#ok<*NASGU>

% handle rectangular patches
if numel(patchSize) == 1
    patchSize = [patchSize patchSize];
end

t0 = tic;
tLast = tic;

Fpower = zeros(patchSize);
count = 0;
for i = 1:length(imageNames)
    if toc(t0) > params.progressStart
        if toc(tLast) > params.progressEvery
            disp([mfilename ' processing image ' int2str(i) '/' int2str(length(imageNames)) ...
                ', ' imageNames(i).path ...
                '... ' num2str(toc(t0), '%.2f') ' seconds elapsed.']);
            tLast = tic;
        end
    end
    
    image = preprocessImage(fullfile(path, imageNames(i).path), blockAvgFactor);
    
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
    
    % XXX the following only works with square patches:
%     imgPatches = patchifyImage(image, patchSize(1));
    
    % Fourier transform the patches
%     imgpatchesF=zeros(size(imgPatches));
%     for j=1:size(imgPatches,3)
%         imgpatchesF(:,:,j)=fft2(imgPatches(:,:,j));
%     end
%     Fpower=Fpower+mean(imgpatchesF.*conj(imgpatchesF),3)/length(images);
end
Fpower = Fpower / count;
Ffilter = 1 ./ sqrt(Fpower);

end