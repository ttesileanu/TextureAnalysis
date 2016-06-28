function [Ffilter] = generateFourierWhitenFilter(imgDirectory, images, blockAvgFactor, patchSize, PCfrac)
% generateFourierWhitenFilter Generate a whitening filter.
%   Ffilter = generateFourierWhitenFilter(imgDirectory, images, blockAvgFactor, ...
%               patchSize, PCfrac)
%   generates a whitening filter using the 'images' from the 'imgDirectory'.
%   The images are first downsampled using 'blockAvgFactor', then split
%   into patches of (linear) size given by 'patchSize'.
%
%   The returned matrix, Ffilter, is the matrix of inverse square roots of
%   the Fourier power values.

% XXX what exactly should PCfrac do?

Fpower=zeros(patchSize);
d=patchSize^2;
numWide = zeros(length(images));
numHigh = zeros(length(images));
for i = 1:length(images)
    LUM_Image = loadLUMImage(fullfile(imgDirectory, images(i).path));
    
    numWide(i)=floor(size(LUM_Image,2)/(blockAvgFactor*patchSize));
    numHigh(i)=floor(size(LUM_Image,1)/(blockAvgFactor*patchSize));
    
    imgpatches=zeros(patchSize,patchSize,numWide(i)*numHigh(i));
    imgpatchesF=zeros(size(imgpatches));
    
    %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
    LUM_Image(LUM_Image <= 0) = 1e-6;
    img=log(LUM_Image);
    
    img = blockAverage(img,blockAvgFactor,'avg');
    % get rid of incomplete patches
    img=img(1:patchSize*numHigh(i),1:patchSize*numWide(i));
    imgres=reshape(img,patchSize,numHigh(i),patchSize*numWide(i));
    for j=1:numHigh(i)
        imgpatches(:,:,((j-1)*numWide(i)+1):(j*numWide(i))) = reshape(imgres(:,j,:),patchSize,patchSize,numWide(i));
    end
    if (~(PCfrac==1))
        imgpatchesV=reshape(imgpatches,size(imgpatches,1)*size(imgpatches,2),size(imgpatches,3));
        % XXX what is PCV?
        wts = PCV' * imgpatchesV;
        wts_reduce = wts(floor((1-PCfrac)*d)+1:d,:);
        imgpatchesV = PCV(:,floor((1-PCfrac)*d)+1:d) * wts_reduce;
        imgpatches = reshape(imgpatchesV,size(imgpatches,1),size(imgpatches,2),size(imgpatches,3));
    end
    for j=1:size(imgpatches,3)
        imgpatchesF(:,:,j)=fft2(imgpatches(:,:,j));
    end
    Fpower=Fpower+mean(imgpatchesF.*conj(imgpatchesF),3)/length(images);
end
Ffilter=1./sqrt(Fpower);
end