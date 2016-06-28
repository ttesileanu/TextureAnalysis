function [ covM, ev, imageCoordinates, Ffilter, sharpness ] = analyzeImageSetModNoPC( imgNamesFile, imgDirectory, imgNamesFile_filter, imgDirectory_filter, blockAvgFactor, PCfrac, patchSize)
% analyzeImageSetModNoPC Split images into patches and calculate texture
% statistics per patch.
%   [covM, ev] = analyzeImageSetModNoPC(
%           imgNamesFile, imgDirectory, imgNamesFile_filter, imgDirectory_filter,
%           blockAvgFactor, PCfrac, patchSize)
%   uses the images identified by 'imgNamesFile_filter' in folder
%   'imgDirectory_filter' to generate a whitening filter, then uses this
%   filter to whiten and analyze the images identified by 'imgNamesFile' in
%   folder 'imgDirectory'. The analysis splits the images into patches of 
%   size 'patchSize' after having downsampled them by 'blockAvgFactor', and
%   calculates texture statistics for each patch. Each row in the matrix
%   'ev' contains the values for the 10 independent texture parameters (see
%   getStatistics and processBlock) for each image patch. 'covM' is the
%   covariance matrix for these statistics across the image patches.
%
%   [covM, ev, imageCoordinates, Ffilter, sharpness] = analyzeImageSetModNoPC(...)
%   also returns information regarding the image position for each patch
%   ('imageCoordinates'), the whitening filter ('Ffilter'), and per-patch
%   sharpness information ('sharpness').
%
%   See also: processBlock, getStatistics.

% XXX what should PCfrac do?

images_filter = parseImageNameFile(imgNamesFile_filter, imgDirectory_filter);
images = parseImageNameFile(imgNamesFile, imgDirectory);

tic

%generate filter (1./sqrt(Fourier Power)) using all images (natural and segmented)
disp('Whitening...');
[Ffilter] = generateFourierWhitenFilter(imgDirectory_filter, images_filter, blockAvgFactor, patchSize, PCfrac);

%compute statistics using only segmented images (a subset of those used to build the filter)
disp('Getting Statistics...');
[ev,N,numWide,numHigh, sharpness]=getImgStats(imgDirectory, images, blockAvgFactor, patchSize, PCfrac, Ffilter);

covM=cov(ev);

%save the locations of patches
%disp('Copying image patches...');
for i=1:length(images)
    imageCoordinates.name{i} = images(i).path;
    for j = 1:numWide(i)
        for k = 1:numHigh(i)
            l=(i-1)*N(i)+j+(k-1)*numWide(i); %index to use
            imageCoordinates.image(l) = i;
            imageCoordinates.x(l) = j;
            imageCoordinates.y(l) = k;
        end
    end
end

t = toc;
disp([mfilename ' took ' num2str(t, 2) ' seconds to run.']);

end