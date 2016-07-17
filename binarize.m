function [img_out] = binarize(img,gamma, patchSize)
% binarize Binarize images.
%   img_out = binarize(img, gamma) binarizes 'img' so that the fraction of
%   'on' pixels is (1+gamma)/2.
%
%   img_out = binarize(img, gamma, patchSize) binarizes each pixel by using
%   a threshold calculated from a square patch surrounding the pixel of
%   linear size 'patchSize'. Patches close to the edges are truncated.

if (gamma < -1) || (gamma > 1)
    error([mfilename ':badgamma'], 'gamma must lie between -1 and 1');
end

p = (1+gamma)/2;

if nargin < 3
    threshold = quantile(img(:),p);
    img_out = double(img >= threshold);
else
    img_out = zeros(size(img));
    
    dr = patchSize/2;

    [ymax,xmax]=size(img);
    for i=1:ymax
        for j=1:xmax
            patch = img( max(1,i-dr):min(ymax,i+dr) , max(1,j-dr):min(xmax,j+dr) );
            threshold = quantile(patch(:),p);
            if img(i,j)>threshold
                img_out(i,j)=1;
            end
        end
    end
end

end