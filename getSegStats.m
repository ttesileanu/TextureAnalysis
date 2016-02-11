function [evB,evF] = getSegStats(blockAvgFactor, patchSize, Rblur, blurFG)
      
% A Hermundstad, Dec 2015
% loads and preprocesses images (block-average, whiten, binarize)
% separately estimates statistics in FG / BG
% uses all non-boundary pixels within FG/BG (so FG/BG stats are estimated 
% from different numbers of pixels)

    if nargin<4
        Rblur=0;
        blurFG=0;
    end
    
    if ~ismember(blockAvgFactor,[1,2])
        disp('Error: block average factor outside range')
    end
    
    if ~ismember(patchSize,32)
        disp('Error: patch size outside range')
    end
        
    load('ImgSegs.mat')
    load(sprintf('Ffilter%dx%d.mat',blockAvgFactor,patchSize))
    
    
    imgDirectory='/Users/ann/Research/projects/textures/BerkeleyImageDB/Images/BSDS/';
    imgNamesFile='./BSR_Test_Index.txt';
    
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

    
    
    
    cshift = patchSize/2-1;
    etrim = 10; %trim edges of images
    
    evB=[];evF=[];
    m=1;
    for i = 1:length(images)
        disp(i)
        
        %select good segmentations
        if ~isnan(seg_select(i))
            
            %load image
            load(strcat(imgDirectory,images(i).path)); 
            LUM_Image(LUM_Image<=0) = 10^(-.356); %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
            
            
            %block average image and semgentation matrix
            img = log(LUM_Image);
            img = blockAverage(img,blockAvgFactor,'avg');
            
            m=m+1;
            if Rblur>0 
                img_blur = log(imgaussfilt(LUM_Image,Rblur));
            else
                img_blur = log(LUM_Image);
            end
            
            seg = segs(i).FG(seg_select(i)).fgMat;
            seg = blockAverage(seg,blockAvgFactor,'avg');

            %filter, trim, and binarize image
            img = imfilter(img,circshift(ifft2(Ffilter),[cshift,cshift]));
            img = img(etrim:end-(etrim-1),etrim:end-(etrim-1));
            img_bin = binarize(img,0,patchSize);
            
            img_blur = blockAverage(img_blur,blockAvgFactor,'avg');
            img_blur = imfilter(img_blur,circshift(ifft2(Ffilter),[cshift,cshift]));
            img_blur = img_blur(etrim:end-(etrim-1),etrim:end-(etrim-1));
            img_blur_bin = binarize(img_blur,0,patchSize);

            seg = seg(etrim:end-(etrim-1),etrim:end-(etrim-1));

            %get foreground stats:
            mask=ones(size(seg));
            mask(seg<.5)=NaN;
            if blurFG
                block = mask.*img_blur_bin;
            else
                block = mask.*img_bin;
            end
            [~, ev, ~, ~] = processBlock(block);
            evF = [evF;ev'];

            %get background stats:
            mask=ones(size(seg));
            mask(seg>=.5)=NaN;
            if blurFG
                block = mask.*img_bin;
            else
                block = mask.*img_blur_bin;
            end
            [~, ev, ~, ~] = processBlock(block);
            evB = [evB;ev'];
            
        end
    end
    clear LUM_Image;
end



function [img_out] = binarize(img,gamma,patchSize)
% Binarize img so that the fraction of 'on' pixels is (1+gamma)/2

img_out = zeros(size(img));
if (gamma < -1) || (gamma > 1)
    display('gamma must lie between -1 and 1');
    return;
end

dr = patchSize/2;
p = (1+gamma)/2;
[ymax,xmax]=size(img);
for i=1:size(img,1)
    for j=1:size(img,2)
        patch = img( max(1,i-dr):min(ymax,i+dr) , max(1,j-dr):min(xmax,j+dr) );
        threshold = quantile(patch(:),p);
        if img(i,j)>threshold
            img_out(i,j)=1;
        end
    end
end

end



