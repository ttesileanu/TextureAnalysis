function [evs,N,numWide,numHigh, sharpness] = getImgStats(imgDirectory, images, blockAvgFactor, patchSize, PCfrac, Ffilter)
% getImgStats Calculate image statistics.
%   [evs, N, numWide, numHigh, sharpness] = getImgStats(...
%       imgDirectory, images, blockAvgFactor, patchSize, PCfrac, Ffilter)
%   calculates texture statistics for the 'images' in 'imgDirectory',
%   whitening with 'Ffilter', downsampling by 'blockAvgFactor', and using
%   patches of (linear) size 'patchSize'.
%
%   The function returns:
%    evs:
%       Values of the 10 independent texture parameters for each of the
%       image patches. (see getStatistics and processBlock)
%    N:
%       Vector of number of patches per image.
%    numWide, numHigh:
%       Vectors giving the number of patches that fit horizontally and
%       vertically, respectively, on each image.
%    sharpness:
%       Vector of patch sharpness values. This is calculated using the
%       median of the Laplacian-filtered images.

% XXX what does PCfrac do?

evs=[];
d=patchSize^2;
sharpness = [];

numWide = zeros(length(images));
numHigh = zeros(length(images));
N = zeros(length(images));

for i = 1:length(images)
    % XXX most of the code here is duplicated from generateFourierWhitenFilter
    LUM_Image = loadLUMImage(fullfile(imgDirectory, images(i).path));
    
    numWide(i)=floor(size(LUM_Image,2)/(blockAvgFactor*patchSize));
    numHigh(i)=floor(size(LUM_Image,1)/(blockAvgFactor*patchSize));
    
    N(i)=numWide(i)*numHigh(i); %number of patches per image
    imgpatches=zeros(patchSize,patchSize,N(i));
    
    %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
    % XXX this isn't such a small value; how was it chosen?
    LUM_Image(LUM_Image <= 0) = 10^(-.356);
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
    % also estimating a measure of sharpness based on a Laplacian filter
    lapflt = [0 1 0; 1 -4 1; 0 1 0];
    crt_sharp = zeros(size(imgpatches, 3), 1);
    for j=1:size(imgpatches,3)
        % estimate sharpness for each patch
        filtered = conv2(imgpatches(:, :, j), lapflt, 'valid');
        crt_sharp(j) = median(abs(filtered(:)));
        % binarize
        imgpatches(:,:,j)= binarize(real(ifft2(fft2(imgpatches(:,:,j)) .* Ffilter)),0);
    end
    [~, ev] = getStats(imgpatches);
    evs = [evs;ev]; %#ok<AGROW>
    sharpness = [sharpness; crt_sharp(:)]; %#ok<AGROW>
    %evs((1+(i-1)*N):(i*N),:)=ev;
end

end