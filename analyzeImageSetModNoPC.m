function [ covM, ev, imageCoordinates, Ffilter, sharpness ] = analyzeImageSetModNoPC( imgNamesFile, imgDirectory, imgNamesFile_filter, imgDirectory_filter, blockAvgFactor, PCfrac, patchSize)


%get images for filter
images_filter=struct('path',{});
n=1;
fid=fopen(imgNamesFile_filter);
while (~feof(fid))
    c=fscanf(fid,'%c',1);
    if (c=='=')
        s=fscanf(fid,'%s\n',1); %This may have to be changed for files made in MS notepad. Carriage return is \r\n
        s=strcat(s(1:(length(s)-4)),'_LUM.mat');
        images_filter(n).path=s;
        n=n+1;
    end
end


%get images for estimating statistics
images=struct('path',{});
n=1;
fid=fopen(imgNamesFile);
while (~feof(fid))
    c=fscanf(fid,'%c',1);
    if (c=='=')
        s=fscanf(fid,'%s\n',1); %This may have to be changed for files made in MS notepad. Carriage return is \r\n
        s=strcat(s(1:(length(s)-4)),'_LUM.mat');
        images(n).path=s;
        n=n+1;
    end
end


tic

%generate filter (1./sqrt(Fourier Power)) using all images (natural and segmented)
disp('Whitening...');
[Ffilter] = generateFourierWhitenFilter(imgDirectory_filter, images_filter, blockAvgFactor, patchSize, PCfrac);

%compute statistics using only segmented images (a subset of those used to build the filter)
disp('Getting Statistics...');
[ev,N,numWide,numHigh, sharpness]=getImgStats(imgDirectory, images, blockAvgFactor, patchSize, PCfrac, Ffilter);
covM=cov(ev);

%save the location of patches
disp('Copying image patches...');
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

toc

end
    

    
    function [Ffilter] = generateFourierWhitenFilter(imgDirectory, images, blockAvgFactor, patchSize, PCfrac)
        Fpower=zeros(patchSize);
        d=patchSize^2;
        for i = 1:length(images)
            load(strcat(imgDirectory,images(i).path)); 
            pic1=LUM_Image; %we use images of all the same size
            numWide(i)=floor(size(pic1,2)/(blockAvgFactor*patchSize)); 
            numHigh(i)=floor(size(pic1,1)/(blockAvgFactor*patchSize));
            clear pic1;
            
            imgpatches=zeros(patchSize,patchSize,numWide(i)*numHigh(i));
            imgpatchesF=zeros(size(imgpatches));
        
            ind=find(LUM_Image<=0);
            LUM_Image(ind) = 10^(-6); %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
            img=log(LUM_Image);
            img = blockAverage(img,blockAvgFactor,'avg');
            img=img(1:patchSize*numHigh(i),1:patchSize*numWide(i));
            imgres=reshape(img,patchSize,numHigh(i),patchSize*numWide(i));
            for j=1:numHigh(i)
                imgpatches(:,:,((j-1)*numWide(i)+1):(j*numWide(i))) = reshape(imgres(:,j,:),patchSize,patchSize,numWide(i));
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
    
    function [evs,N,numWide,numHigh, sharpness] = getImgStats(imgDirectory, images, blockAvgFactor, patchSize, PCfrac, Ffilter)
        
        evs=[];
        d=patchSize^2;
        sharpness = [];
        for i = 1:length(images)
            load(strcat(imgDirectory,images(i).path)); 
            
            pic1=LUM_Image; %we use images of all the same size
            numWide(i)=floor(size(pic1,2)/(blockAvgFactor*patchSize)); 
            numHigh(i)=floor(size(pic1,1)/(blockAvgFactor*patchSize));
            clear pic1;
            
            N(i)=numWide(i).*numHigh(i); %number of patches per image
            imgpatches=zeros(patchSize,patchSize,N(i));
            
            ind=find(LUM_Image<=0);
            LUM_Image(ind) = 10^(-.356); %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
            img=log(LUM_Image);
            img = blockAverage(img,blockAvgFactor,'avg');
            img=img(1:patchSize*numHigh(i),1:patchSize*numWide(i));
            imgres=reshape(img,patchSize,numHigh(i),patchSize*numWide(i));
            for j=1:numHigh(i)
                imgpatches(:,:,((j-1)*numWide(i)+1):(j*numWide(i))) = reshape(imgres(:,j,:),patchSize,patchSize,numWide(i));
            end
            if (~(PCfrac==1))
                imgpatchesV=reshape(imgpatches,size(imgpatches,1)*size(imgpatches,2),size(imgpatches,3));
                wts = PCV' * imgpatchesV;
                wts_reduce = wts(floor((1-PCfrac)*d)+1:d,:);
                imgpatchesV = PCV(:,floor((1-PCfrac)*d)+1:d) * wts_reduce;
                imgpatches = reshape(imgpatchesV,size(imgpatches,1),size(imgpatches,2),size(imgpatches,3));
            end
            lapflt = [0 1 0; 1 -4 1; 0 1 0];
            crt_sharp = zeros(size(imgpatches, 3), 1);
            for j=1:size(imgpatches,3)
                % estimate sharpness for each patch
                filtered = conv2(imgpatches(:, :, j), lapflt, 'valid');
                crt_sharp(j) = median(abs(filtered(:)));
                % binarize
                imgpatches(:,:,j)= binarize(real(ifft2(fft2(imgpatches(:,:,j)) .* Ffilter)),0);
            end
            [P ev entropy] = getStats(imgpatches);
            evs = [evs;ev];
            sharpness = [sharpness; crt_sharp(:)];
            %evs((1+(i-1)*N):(i*N),:)=ev;
            clear LUM_Image;
        end
    end