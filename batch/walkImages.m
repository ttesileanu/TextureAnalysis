function results = walkImages(pipeline, images, varargin)
% walkImages Apply a series of functions to every image in a set and
% collect the results.
%   results = walkImages(pipeline, images) loads every image identified
%   by the `images` cell array and passes it through a series of functions
%   from the cell array `pipeline`. It collects the results from the last
%   function application into a structure.
%
%   Each element of the `pipeline` except the last should be a function
%   that takes a matrix or a 3d array of patches and returns two arguments:
%   another image (2d or 3d), and a matrix of cropping/scaling information
%   containing rows of the form [row1, col1, row2, col2]:
%       [imgOut, crop] = fct(img)
%   Two kinds of functions are allowed
%       (a) functions that act on each patch, potentially removing some. In
%           this case, the number of patches in the output is at most as
%           large as that in the input (i.e., `size(imgOut, 3) <= size(img, 3)`),
%           and each row in `crop` matches one patch in `img` (i.e.,
%           `size(crop, 1) == size(img, 3)`). Rows that are all equal to 0
%           indicate patches that have been dropped in the processing. The
%           number of such rows must thus be equal to the difference
%           between the number of patches in the input and output images.
%           Each row indicates a potential crop within each patch.
%       (b) functions that expand a 2d matrix into a 3d array of patches.
%           In this case, the number of rows in `crop` matches the number
%           of patches in `imgOut`, `size(crop, 1) == size(imgOut, 3)`, and
%           each row indicates the extents of each patch in the input image.
%   If the image that is returned is an empty matrix for any function
%   within the `pipeline`, the current image is ignored. Empty crop
%   information implies no cropping compared to the input.
%
%   The last element of the processing `pipeline` should be a function of
%   the form
%       resStruct = fct(i, image, crop)
%   taking the index in the image set of the current image being processed,
%   `i`, the image data itself, `image`, and an nPatches x 4 matrix giving
%   the locations of all the patches in the initial image (where `nPatches`
%   is equal to `size(image, 3)`). The function should return either an
%   empty matrix (in which case the result is discarded), or a structure.
%   All the fields in the structures returned for all the images are
%   combined into the output `results`.
%
%   The `images` argument to walkImages can be input as
%     * file names, in which case they are loaded using `loadLUMImage`;
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
    
    % initial crop -- a single patch covering the whole image
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
        oldNPatches = size(crtImage, 3);
        [crtImage, crtCrop] = pipeline{j}(crtImage);
        newNPatches = size(crtImage, 3);
        
        % an empty result leads to skipping
        if isempty(crtImage)
            skipImage = true;
            break;
        end
        
        % keep old crop if new one is empty
        if ~isempty(crtCrop)
            % new crop -- one row for each patch
            newCrop = zeros(newNPatches, 4);

            if newNPatches > oldNPatches
                % we're patchifying a 2d image
                if oldNPatches ~= 1
                    error([mfilename ':badpatch'], ...
                        'Can only patchify 2d image. (image %d, stage %d)', i, j);
                end
            else
                % we might be losing some of the patches
                keepMask = any(crtCrop ~= 0, 2);
                % select only the patches we're keeping
                crop = crop(keepMask, :);
                crtCrop = crtCrop(keepMask, :);
                % check match between number of patches and rows of crtCrop
                if size(crtCrop, 1) ~= newNPatches
                    error([mfilename ':badnp'], ...
                        'Number of patches in image output doesn''t match non-zero rows of crop.  (image %d, stage %d)', ...
                        i, j);
                end
                % note that newNPatches can't be zero (we checked for empty
                % crtImage earlier), and so crtCrop and crop can also not
                % be empty here
            end
            
            % figure out any scaling between the image coordinates for the
            % input patches at this particular pipeline position, and the
            % coordinates in the initial image
            scaleRows = (crop(:, 3) - crop(:, 1) + 1) ./ oldSize(:, 1);
            scaleCols = (crop(:, 4) - crop(:, 2) + 1) ./ oldSize(:, 2);
            
            % find new patch extents in initial-image coordinates
            newCrop(:, 1) = crop(:, 1) + scaleRows.*(crtCrop(:, 1) - 1);
            newCrop(:, 2) = crop(:, 2) + scaleCols.*(crtCrop(:, 2) - 1);
            newCrop(:, 3) = newCrop(:, 1) + scaleRows.*(crtCrop(:, 3) - crtCrop(:, 1) + 1) - 1;
            newCrop(:, 4) = newCrop(:, 2) + scaleCols.*(crtCrop(:, 4) - crtCrop(:, 2) + 1) - 1;

            % update the crop
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
    else
        if ~params.quiet
            if ~generate && ischar(images{i})
                crtName = ['(' images{i} ')'];
            else
                crtName = '';
            end
            warning('Image #%d %s skipped because of empty output at stage %d.', ...
                i, crtName, j);
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