function [Ffilter] = generateFourierWhitenFilter_OA(imgDirectory, images, blockAvgFactor, patchSize, PCfrac, numWide, numHigh)
    Fpower=zeros(patchSize);
    imgpatches=zeros(patchSize,patchSize,numWide*numHigh);
    imgpatchesF=zeros(size(imgpatches));
    d=patchSize^2;
    for i = 1:length(images)
        load(strcat(imgDirectory,images(i).path)); 
        ind=find(LUM_Image<=0);
        LUM_Image(ind) = 10^(-6); %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
        img=log(LUM_Image);
        img = blockAverage(img,blockAvgFactor,'avg');
        img=img(1:patchSize*numHigh,1:patchSize*numWide);
        imgres=reshape(img,patchSize,numHigh,patchSize*numWide);
        for j=1:numHigh
            imgpatches(:,:,((j-1)*numWide+1):(j*numWide)) = reshape(imgres(:,j,:),patchSize,patchSize,numWide);
        end
        if (~(PCfrac==1))
            imgpatchesV=reshape(imgpatches,size(imgpatches,1)*size(imgpatches,2),size(imgpatches,3));
            wts = PCV' * imgpatchesV;
            wts_reduce = wts(floor((1-PCfrac)*d)+1:d,:);
            imgpatchesV = PCV(:,floor((1-PCfrac)*d)+1:d) * wts_reduce;
            imgpatches = reshape(imgpatchesV,size(imgpatches,1),size(imgpatches,2),size(imgpatches,3));
        end
        for j=1:size(imgpatches,3)
            imgpatchesF(:,:,j)=fft2(imgpatches(:,:,j));
        end
        Fpower=Fpower+mean(imgpatchesF.*conj(imgpatchesF),3)/length(images);

        clear LUM_Image;
    end
    Ffilter=1./sqrt(Fpower);
end