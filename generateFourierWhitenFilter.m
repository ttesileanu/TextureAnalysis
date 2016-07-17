function [Ffilter] = generateFourierWhitenFilter(imgDirectory, images, blockAvgFactor, patchSize)
% generateFourierWhitenFilter Generate a whitening filter.
%   Ffilter = generateFourierWhitenFilter(imgDirectory, images, ...
%               blockAvgFactor, patchSize)
%   generates a whitening filter using the 'images' from the 'imgDirectory'.
%   The images are first downsampled using 'blockAvgFactor', then split
%   into patches of (linear) size given by 'patchSize'.
%
%   The returned matrix, Ffilter, is the matrix of inverse square roots of
%   the Fourier power values.

Fpower=zeros(patchSize);
for i = 1:length(images)
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