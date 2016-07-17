function [evs,numWide,numHigh, others] = ...
    getImgStats(...
        imgDirectory, imageNames, ...
        blockAvgFactor, patchSize, Ffilter, varargin)
% getImgStats Calculate image statistics.
%   [evs, N, numWide, numHigh, imageCoordinates, sharpness] = getImgStats(...
%       imgDirectory, imageNames, blockAvgFactor, patchSize, Ffilter)
%   calculates texture statistics for the images with the given 'imageNames'
%   found in 'imgDirectory', downsampling by 'blockAvgFactor', whitening
%   with 'Ffilter', and using patches of (linear) size 'patchSize'.
%
%   getImgStats(..., 'segmentation', segmentation) uses the segmentation
%   data from the 'segmentation' structure to focus on particular areas of
%   the images. Note that the analysis will still return the full
%   statistics in 'evs', but will return statistics for areas within the
%   segmentation masks in 'others.evF' and for areas outside the
%   segmentation masks in 'others.evB'. The 'segmentation' structure needs
%   to have fields 'segs', 'segSel', 'segField', and 'segFct; for a
%   description of the meaning of these fields, see processData.m. Images
%   for which 'segSel' is NaN are ignored.
%
%   The function returns:
%    evs:
%       Values of the 10 independent texture parameters for each of the
%       image patches. (see getStatistics and processBlock)
%    numWide, numHigh:
%       Vectors giving the number of patches that fit horizontally and
%       vertically, respectively, on each image.
%    others:
%       A structure containing other detailed information regarding the
%       analysis. This contains the following fields:
%        'imageCoordinates':
%           A structure containing information about the location and image
%           of origin for each of the patches that were analyzed. This is
%           provided in fields 'image', 'x', and 'y'. It also contains a
%           copy of the image names, images(i).path, in a cell array called
%           'name'.
%        'sharpness':
%           Vector of patch sharpness values. This is calculated using the
%           median of the Laplacian-filtered images.
%        'origImages', 'blockAFImages':
%           Cell arrays containing the image matrices for the original
%           grayscale image, and the block-averaged image, respectively.
%           These are present only if the 'images' option is true.
%        'whitenedImages', 'binarizedImages':
%           Cell arrays containing the images after the whitening filter,
%           and after whitening and binarization, respectively. These are
%           present only if the 'images' option is true.
%
%   Options:
%    'images': bool
%       If true, the details structure contains outputs 'origImages',
%       'blockAFImages', 'whitenedImages', and 'binarizedImages' (see
%       above).
%    'filterFull': bool
%       If true, the whitening filter is applied to the full image. If
%       false, it is applied to each patch separately.
%    'fullImageEv': bool
%       If true, calculate a single set of 10 'ev' values for the whole
%       image (and another two sets 'evF' and 'evB' if a segmentation is
%       present). For now this only works when 'filterFull' is true.
%    'overlappingPatches': bool
%       Set to true to use a pixel-by-pixel sliding patch to evaluate
%       statistics. This generates many more patches but these are no
%       longer independent of each other.
%    'segmentation': struct
%       If this is provided, then the analysis will not only be run on the
%       full images in the dataset, but also on those subsets of the images
%       that are inside the masks provided by the segmentation (see
%       processData.m), and also for the complementary masks. The values
%       will be returned as 'others.evF', 'others.covF', and 'others.evB',
%       'others.covB', respectively.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('segmentation', [], @(s) isempty(s) || isstruct(s));
parser.addParameter('images', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('filterFull', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('fullImageEv', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('overlappingPatches', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results; %#ok<*NASGU>

segmentation = params.segmentation;

evs=[];
sharpness = [];

numWide = zeros(length(imageNames), 1);
numHigh = zeros(length(imageNames), 1);

imageCoordinates.name = cell(1, length(imageNames));
imageCoordinates.image = [];
imageCoordinates.x = [];
imageCoordinates.y = [];

others = struct;
if params.images
    others.origImages = cell(1, length(imageNames));
    others.blockAFImages = cell(1, length(imageNames));
    others.whitenedImages = cell(1, length(imageNames));
    others.binarizedImages = cell(1, length(imageNames));
end

if ~isempty(segmentation) && ~isempty(segmentation.segs)
    % create segmentation-based outputs
    others.evF = [];
    others.evB = [];
end

patchSizeOrig = patchSize*blockAvgFactor;
if params.filterFull
    % will need to trim if we do the full-image filter
    trimEdge = patchSize/2;
    %trimEdge = 10;
else
    trimEdge = 0;
end
trimEdgeOrig = trimEdge*blockAvgFactor;

for i = 1:length(imageNames)
    imageCoordinates.name{i} = imageNames(i).path;
    
    [preprocessed, image0] = preprocessImage(fullfile(imgDirectory, imageNames(i).path), ...
        blockAvgFactor);
    
    % are we to apply the whitening filter to the full image?
    if params.filterFull
        % yes, filter and binarize on full image
        cshift = patchSize/2-1;
        whitenedImage = imfilter(preprocessed, circshift(ifft2(Ffilter), [cshift,cshift]));
        
        whitenedImage = whitenedImage(trimEdge:end-trimEdge+1, trimEdge:end-trimEdge+1);
                
        % XXX this is SLOW!
        binarizedImage = binarize(whitenedImage, 0, patchSize);
        
        if ~params.overlappingPatches && ~params.fullImageEv
            [imgPatches, numWide(i), numHigh(i), xcoords, ycoords] = ...
                patchifyImage(binarizedImage, patchSize);
            finalWidth = numWide(i)*patchSize;
            finalHeight = numHigh(i)*patchSize;
            
            whitenedImage = whitenedImage(1:finalHeight, 1:finalWidth);
            binarizedImage = binarizedImage(1:finalHeight, 1:finalWidth);
        else
            finalHeight = size(binarizedImage, 1);
            finalWidth = size(binarizedImage, 2);
        end
    else
        % no, filter and binarize per-patch
        [imgPatches, numWide(i), numHigh(i), xcoords, ycoords] = ...
            patchifyImage(preprocessed, patchSize);
        finalWidth = numWide(i)*patchSize;
        finalHeight = numHigh(i)*patchSize;

        binarizedImage = zeros(finalHeight, finalWidth);
        whitenedImage = zeros(finalHeight, finalWidth);
        for j=1:size(imgPatches,3)
            xr = (xcoords(j) - 1)*patchSize + 1;
            yr = (ycoords(j) - 1)*patchSize + 1;
            
            whitenedPatch = real(ifft2(fft2(imgPatches(:,:,j)) .* Ffilter));
            whitenedImage(yr:yr+patchSize-1, xr:xr+patchSize-1) = whitenedPatch;
            
            binarizedPatch = binarize(whitenedPatch,0);
            imgPatches(:,:,j) = binarizedPatch;
            binarizedImage(yr:yr+patchSize-1, xr:xr+patchSize-1) = binarizedPatch;
        end
    end

    % store images, if request
    if params.images
        others.origImages{i} = image0;
        others.blockAFImages{i} = preprocessed;
        others.whitenedImages{i} = whitenedImage;
        others.binarizedImages{i} = binarizedImage;
    end

    % prepare segmentation, if available
    if ~isempty(segmentation) && ~isempty(segmentation.segs)
        if ~isempty(segmentation.segSel) && isnan(segmentation.segSel(i))
            % skip this image
            continue;
        end
        
        crtSeg = segmentation.segs(i);
        if isfield(crtSeg, 'FG')
            crtSeg = crtSeg.FG;
        end
        if ~isempty(segmentation.segSel)
            crtSeg = crtSeg(segmentation.segSel(i));
        end
        
        crtImageMask = crtSeg.(segmentation.segField);
        if ~isempty(segmentation.segFct)
            crtImageMask = segmentation.segFct(crtImageMask);
        end
        
        % if we used full-image filtering, we had to trim some edges
        if trimEdgeOrig > 0
            crtImageMask = crtImageMask(trimEdgeOrig:end - trimEdgeOrig + 1, ...
                trimEdgeOrig:end - trimEdgeOrig + 1);
        end
    end

    
    if ~params.overlappingPatches && ~params.fullImageEv
        % calculate sharpness per patch
        % Laplacian filter
        lapflt = [0 1 0; 1 -4 1; 0 1 0];
        crtSharp = zeros(size(imgPatches, 3), 1);
        
        if trimEdgeOrig > 0
            image0Trimmed = image0(trimEdgeOrig:end-trimEdgeOrig+1, trimEdgeOrig:end-trimEdgeOrig+1);
        else
            image0Trimmed = image0;
        end
        for j=1:size(imgPatches,3)
            xr = (xcoords(j) - 1)*patchSizeOrig + 1;
            yr = (ycoords(j) - 1)*patchSizeOrig + 1;
            crtPatch = image0Trimmed(yr:yr+patchSizeOrig-1, xr:xr+patchSizeOrig-1);
            
            % estimate sharpness for each patch
            filtered = conv2(crtPatch, lapflt, 'valid');
            crtSharp(j) = median(abs(filtered(:)));
        end
        
        % store origin images and patch positions
        N = numWide(i)*numHigh(i);
        imageCoordinates.image = [imageCoordinates.image ; repmat(i, N, 1)];
        imageCoordinates.x = [imageCoordinates.x ; xcoords(:)];
        imageCoordinates.y = [imageCoordinates.y ; ycoords(:)];
        
        % get and store statistics
        [~, ev] = getStats(imgPatches);
        evs = [evs;ev]; %#ok<AGROW>
        sharpness = [sharpness; crtSharp(:)]; %#ok<AGROW>

        % handle segmentation, if available
        if ~isempty(segmentation) && ~isempty(segmentation.segs)
            % find patches completely contained in the mask, and those
            % completely outside
            fgPatchMask = false(size(xcoords));
            bgPatchMask = false(size(xcoords));
            for l = 1:length(xcoords)
                x0 = (xcoords(l) - 1)*patchSizeOrig + 1;
                y0 = (ycoords(l) - 1)*patchSizeOrig + 1;
                maskPatch = crtImageMask(y0:y0 + patchSizeOrig-1, x0:x0 + patchSizeOrig-1);
                fgPatchMask(l) = all(maskPatch(:));
                bgPatchMask(l) = ~any(maskPatch(:));
            end
            
            fgPatches = imgPatches(:, :, fgPatchMask);
            bgPatches = imgPatches(:, :, bgPatchMask);
            
            if ~isempty(fgPatches)
                [~, evF] = getStats(fgPatches);
                others.evF = [others.evF ; evF];
            end
            if ~isempty(bgPatches)
                [~, evB] = getStats(bgPatches);
                others.evB = [others.evB ; evB];
            end
            
            % XXX should calculate&store imageCoordinates for segmented images,
            % too
        end
    elseif params.overlappingPatches
        % overlappingPatches is true
        nCols = finalWidth - patchSize + 1;
        nRows = finalHeight - patchSize + 1;
        nPatches = nCols*nRows;
        crtEv = zeros(nPatches, 10);
        crtEvF = zeros(nPatches, 10);
        crtEvB = zeros(nPatches, 10);
        crt = 1;
        crtF = 1;
        crtB = 1;
        
        for y = 1:nRows
            for x = 1:nCols
                crtPatch = binarizedImage(y:y+patchSize-1, x:x+patchSize-1);
                
                [~, ev] = processBlock(crtPatch);
                crtEv(crt, :) = ev;
                if ~isempty(segmentation) && ~isempty(segmentation.segs)
                    % mask coordinates can be different
                    xOrig = 1 + (x-1)*blockAvgFactor;
                    yOrig = 1 + (y-1)*blockAvgFactor;
                    crtMask = crtImageMask(yOrig:yOrig+patchSizeOrig-1, ...
                        xOrig:xOrig+patchSizeOrig-1);
                    if all(crtMask(:))
                        crtEvF(crtF, :) = ev;
                        crtF = crtF + 1;
                    elseif ~any(crtMask(:))
                        crtEvB(crtB, :) = ev;
                        crtB = crtB + 1;
                    end
                end
                
                crt = crt + 1;
            end
        end
        evs = [evs ; crtEv]; %#ok<AGROW>
        crtEvF = crtEvF(1:crtF-1, :);
        crtEvB = crtEvB(1:crtB-1, :);
        others.evF = [others.evF ; crtEvF];
        others.evB = [others.evB ; crtEvB];
    elseif params.fullImageEv
        % treat the whole image as a single block
        [~, ev] = processBlock(binarizedImage);
        evs = [evs; ev(:)']; %#ok<AGROW>
        
        smallMask = blockAverage(crtImageMask, blockAvgFactor, 'wta');
        imgF = binarizedImage;
        % using NaNs to ignore pixels outside of mask
        imgF(smallMask) = nan;
        [~, evF] = processBlock(imgF);
        others.evF = [others.evF ; evF(:)'];
        
        imgB = binarizedImage;
        imgB(~smallMask) = nan;
        [~, evB] = processBlock(imgB);
        others.evB = [others.evB ; evB(:)'];
    end
end

others.imageCoordinates = imageCoordinates;
others.sharpness = sharpness;

end