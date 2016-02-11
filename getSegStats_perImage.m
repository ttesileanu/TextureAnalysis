function [evB,evF] = getSegStats_perImage(nImg, R, blockAvgFactor, patchSize, Rblur, blurFG)


    % A Hermundstad, Feb 2016
    % Loads and preprocesses single image (block-average, whiten,
    % binarize).  Uses sliding square glider to separately estimate 
    % statistics in FG/BG patches (note that there are unequal numbers 
    % of FG vs BG patches)

    %nImg = index of image to use
    %R = size of sliding quare glider (typically R=16,24,32,...)
    if nargin<6
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
    
    evB=[];evF=[];
    
    cshift = patchSize/2-1;
    etrim = 10;
    for i = nImg

        if ~isnan(seg_select(i))
            
            %load image
            load(strcat(imgDirectory,images(i).path)); 

            LUM_Image(LUM_Image<=0) = 10^(-.356); %Set any luminance values which have been corrected to below zero to a very small number so we can take the logarithm
            img = log(LUM_Image);
            
            if Rblur>0 
                img_blur = log(imgaussfilt(LUM_Image,Rblur));
            else
                img_blur = log(LUM_Image);
            end
            
            
            img = blockAverage(img,blockAvgFactor,'avg');
            img_blur = blockAverage(img_blur,blockAvgFactor,'avg');

            seg = segs(i).FG(seg_select(i)).fgMat;
            seg = blockAverage(seg,blockAvgFactor,'avg');
            
            img = imfilter(img,circshift(ifft2(Ffilter),[cshift,cshift]));
            img = img(etrim:end-(etrim-1),etrim:end-(etrim-1));
            img_bin = binarize(img,0,patchSize);
            
            img_blur = imfilter(img_blur,circshift(ifft2(Ffilter),[cshift,cshift]));
            img_blur = img_blur(etrim:end-(etrim-1),etrim:end-(etrim-1));
            img_blur_bin = binarize(img_blur,0,patchSize);

            imgF=zeros(size(img));
            
            seg = seg(etrim:end-(etrim-1),etrim:end-(etrim-1));

            %R = patchSize; %4x4 patches
            [lh,lw] = size(img);
            nrows = lh-(R-1);
            ncols = lw-(R-1);
            
            for x=1:nrows
                for y=1:ncols
                    block = img_bin(x:(x+R-1),y:(y+R-1));
                    block_blur = img_blur_bin(x:(x+R-1),y:(y+R-1));
                    m = mean(mean(seg(x:(x+R-1),y:(y+R-1)) ));

                    if m==1 %foreground
                        if blurFG
                            [~, ev, ~, ~] = processBlock(block_blur);
                        else
                            [~, ev, ~, ~] = processBlock(block);
                        end
                        evF = [evF;ev'];

                    elseif m==0 %background
                        if blurFG
                            [~, ev, ~, ~] = processBlock(block);
                        else
                            [~, ev, ~, ~] = processBlock(block_blur);
                        end
                        evB = [evB;ev'];
                    end
                end
            end
            
            

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

function [img_out] = grayscale_convert(img,gamma,patchSize)
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
        threshold = quantile(patch(:),2);
        if img(i,j)>threshold(2)
            img_out(i,j)=1;
        elseif img(i,j)>threshold(1) && img(i,j)<=threshold(2)
            img_out(i,j)=0;
        else
            img_out(i,j)=-1;
        end
    end
end
end



