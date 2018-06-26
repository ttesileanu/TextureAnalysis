function results = walkImages(pipeline, images, varargin)
% walkImages Apply a series of functions to every image in a set and
% collect the results.
%   results = walkImages(pipeline, images) loads every image identified
%   by the `images` cell array and passes it through a series of functions
%   from the cell array `pipeline`. It collects the results from the last
%   function application into a structure.
%
%   Each element of the `pipeline` except the last should be a function
%   that takes one 2d or 3d image (the third dimenion represents patches in
%   the case in which the image was already split into patches) and returns
%   two arguments:
%     (a) in the 2d case, an image, and a cropping/scaling matrix
%         [row1, col1, row2, col2] identifying the region in the input
%         image (for this particular function) corresponding to the
%         resulting image.
%     (b) in the 3d case, a 3d (i.e., patchified) image, and either an
%         npatches x 4 matrix of patch coordinates in the format
%         [row1, col1, row2, col2] (if this is the first function that
%         patchified the image), or an empty matrix to keep the same
%         patches as in the input.
%   If the image that is returned is an empty matrix for any function
%   within the `pipeline`, the current image is ignored. An empty crop or
%   list of patch locations implies using the patches from the previous
%   pipeline function.
%
%   The last element of the processing `pipeline` should be a function of
%   the form
%       fct(i, image, crop)
%   taking the index in the image set of the current image being processed,
%   `i`, the image data itself, `image`, and, when the image is 2d, a
%   cropping/scaling matrix in the format [row1, col1, row2, col2]
%   identifying the region in the initial image (before any of the functions
%   from `pipeline` were applied) corresponding to the final image (the one
%   that's passed to `fct`). If the image is 3d (i.e., patchified), the
%   `crop` argument is an nPatches x 4 matrix giving the locations of all
%   the patches in the initial image. The function should return either an
%   empty matrix (in which case the result is discarded), or a structure.
%   All the fields in the structures returned for each image are combined
%   into the output `results`.
%
%   The `images` can be input as
%     * file names, in which case they are loaded using loadLUMImage;
%     * arrays, in which case they are used as-is;
%     * a pair, `{n, imageGenerator}`, in which case `n` images are
%       generated on the fly by calling the `imageGenerator`; the ith
%       image is generated using `imageGenerator(i)`.
%
%   Options:
%    'quiet'
%       Set to `true` to avoid displaying a progress bar.
%    'storeImages'
%       Set to `true` to keep track of all the images generated while
%       walking the image set, including the output from every pipeline
%       function except the last.
%
%   See also: loadLUMImage.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('quiet', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('storeImages', false, @(b) islogical(b) && isscalar(b));

% display defaults
if nargin == 1 && strcmp(pipeline, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle the various ways of inputting images
if length(images) == 2 && isscalar(images{1}) && isnumeric(images{1}) && ...
        isa(images{2}, 'function_handle')
    generate = true;
    imageCount = images{1};
    imageGenerator = images{2};
else
    generate = false;
    imageCount = length(images);
end

% display a progress bar if allowed
if ~params.quiet
    progress = TextProgress(['image 0/' int2str(imageCount)]);
end

results = [];
for i = 1:imageCount
    % load/generate the image
    if generate
        crtImage = imageGenerator(i);
    else
        if ischar(images{i})
            crtImage = loadLUMImage(images{i});
        else
            crtImage = images{i};
        end
    end
    
    % initial crop covers the whole image
    crop = [1 1 size(crtImage)];
    
    % do the preprocessing (pipeline(1:end-1))
    crtImageCopies = struct;
    skipImage = false;
    for j = 1:length(pipeline)-1
        % store image&crop if asked to
        if params.storeImages
            crtImageCopies.(['stage' int2str(j-1)]) = crtImage;
            crtImageCopies.(['crop' int2str(j-1)]) = crop;
        end
        
        % process
        oldSize = size(crtImage);
        [crtImage, crtCrop] = pipeline{j}(crtImage);
        
        % an empty result leads to skipping
        if isempty(crtImage)
            skipImage = true;
            if ~params.quiet
                if ~generate && ischar(images{i})
                    crtName = ['(' images{i} ')'];
                else
                    crtName = '';
                end
                warning('Image  #%d %s skipped because of empty output at stage #d.', ...
                    i, crtName, j);
            end
            break;
        end
        
        % keep old crop if new one is empty
        if ~isempty(crtCrop)
            % can only generate a new crop if the previous image was 2d
            if ~isvector(crop)
                error([mfilename ':badcrop'], 'Cannot recrop patchified images.');
            end
            % for 2d images, make sure crop is row vector
            if isvector(crtCrop)
                crtCrop = crtCrop(:)';
            end
            % update the crop
            scaleRows = (crop(3) - crop(1) + 1) / oldSize(1);
            scaleCols = (crop(4) - crop(2) + 1) / oldSize(2);
            % I *think* this is right...
            newCrop = zeros(size(crtCrop, 1), length(crop)); 
            newCrop(:, 1) = crop(1) + scaleRows*(crtCrop(:, 1) - 1);
            newCrop(:, 2) = crop(2) + scaleCols*(crtCrop(:, 2) - 1);
            newCrop(:, 3) = crop(1) + scaleRows*(crtCrop(:, 3) - crtCrop(:, 1) + 1) - 1;
            newCrop(:, 4) = crop(2) + scaleCols*(crtCrop(:, 4) - crtCrop(:, 2) + 1) - 1;
            crop = newCrop;
        end
    end
    
    % only continue if we haven't gotten any empty results
    if ~skipImage
        % store last image&crop, if requested
        if params.storeImages
            crtImageCopies.(['stage' int2str(length(pipeline)-1)]) = crtImage;
            crtImageCopies.(['crop' int2str(length(pipeline)-1)]) = crop;
        end
        
        % call the last function, pipeline(end)
        crtResult = pipeline{end}(i, crtImage, crop);
        
        % skip entries for which the result is empty
        if ~isempty(crtResult)
            % add image copies, if requested
            if params.storeImages
                crtResult.imageCopies = crtImageCopies;
            end
            
            % add the current stats to the overall structure
            crtFields = fieldnames(crtResult);
            for j = 1:numel(crtFields)
                crtField = crtFields{j};
                crtValue = crtResult.(crtField);
                % if field already exists in results, concatenate
                if isfield(results, crtField)
                    results.(crtField) = [results.(crtField) ; crtValue];
                else
                    % ...otherwise create the field, making sure to build a
                    % cell array for strings
                    if ischar(crtValue)
                        results.(crtField) = {crtValue};
                    else
                        results.(crtField) = crtValue;
                    end
                end
            end
        end
    end
    
    % update progress bar if allowed
    if ~params.quiet
        progress.update(100*i/imageCount, 'prefix', ['image ' int2str(i) ...
            '/' int2str(imageCount)]);
    end
end

% update progress bar if allowed
if ~params.quiet
    progress.done('done!');
end

if ~generate
    % keep a copy of the image names that were used
    results.imageNames = images;
end
results.imageCount = imageCount;

end