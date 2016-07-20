function [Ffilter] = generateFourierWhitenFilter(imgDirectory, images, blockAvgFactor, patchSize, ...
    varargin)
% generateFourierWhitenFilter Generate a whitening filter.
%   Ffilter = generateFourierWhitenFilter(imgDirectory, images, ...
%               blockAvgFactor, patchSize)
%   generates a whitening filter using the 'images' from the 'imgDirectory'.
%   The images are first downsampled using 'blockAvgFactor', then split
%   into patches of (linear) size given by 'patchSize'.
%
%   The returned matrix, Ffilter, is the matrix of inverse square roots of
%   the Fourier power values.
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

t0 = tic;
tLast = tic;

Fpower=zeros(patchSize);
for i = 1:length(images)
    if toc(t0) > params.progressStart
        if toc(tLast) > params.progressEvery
            disp([mfilename ' processing image ' int2str(i) '/' int2str(length(images)) ...
                ', ' images(i).path ...
                '... ' num2str(toc(t0), '%.2f') ' seconds elapsed.']);
            tLast = tic;
        end
    end
    
    image = preprocessImage(fullfile(imgDirectory, images(i).path), blockAvgFactor);
    imgPatches = patchifyImage(image, patchSize);
    
    % Fourier transform the patches
    imgpatchesF=zeros(size(imgPatches));
    for j=1:size(imgPatches,3)
        imgpatchesF(:,:,j)=fft2(imgPatches(:,:,j));
    end
    Fpower=Fpower+mean(imgpatchesF.*conj(imgpatchesF),3)/length(images);
end
Ffilter=1./sqrt(Fpower);
end