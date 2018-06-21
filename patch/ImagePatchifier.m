classdef ImagePatchifier < handle
    % ImagePatchifier Traverse an image and generate patches.
    %
    % The function can be used to generate equally-spaced square or
    % rectangular patches that are either overlapping or non-overlapping.
    %
    % Note: the function keeps a copy of the image. Because of Matlab's
    % lazy-copying policy, this does not take up memory *unless* you edit
    % the image while the patchifier is in use.
    %
    % Examples:
    %
    %   1. Non-overlapping patches
    %
    %   patchifier = ImagePatchifier(image, 32);
    %
    %   means = zeros(floor(size(image)/32));
    %   while patchifier.next
    %       patch = patchifier.get;
    %       means(patchifier.gridIndex(1), patchifier.gridIndex(2)) = ...
    %           mean(patch(:));
    %   end
    %
    %   2. Overlapping, rectangular patches
    %
    %   image = randn(24);
    %   patchifier = ImagePatchifier(image, [8 16], [3 5]);
    %   disp(['There are ' int2str(patchifier.numPatches) ' patches.']);
    %   disp('Here are their locations in the original image:');
    %   disp(patchifier.patchLocations);
    %
    %   3. Need to call 'next' before grabbing patches
    %
    %   patchifier = ImagePatchifier(randn(24), 8);
    %   try
    %       imagesc(patchifier.get);
    %   catch
    %       disp('Can''t get patch before calling ''next''. ');
    %   end
    %
    %   4. Current patch region
    %
    %   patchifier = ImagePatchifier(randn(32), 16);
    %   patchifier.next;
    %   disp('First patch at ');
    %   disp(patchifier.getPatchCoordinates);
    %
    %   5. Generate coordinates only
    %
    %   patchifier = ImagePatchifier({[33 17]}, 8);
    %   disp(patchifier.patchLocations);
    
    properties
        image = [];         % image we're working on
        imSize = [];        % size of image we're working on
        patchSize = [];     % patch size we're using
        stride = [];        % stride we're using
        gridSize = [];      % number of patches in each direction
%         patchLocations = [];% numPatches x 2 matrix of patch locations in image
%         gridIndices = [];   % numPatches x 2 matrix of grid indices
        gridIndex = [];     % row&column index of current grid position
    end
    
    properties (Access=private)
        current_ = [];      % current index in patchLocations vector
    end
    
    methods
        function obj = ImagePatchifier(image, patchSize, stride)
            % ImagePatchifier Constructor.
            %   ImagePatchifier(image, patchSize) uses non-overlapping
            %   patches of the given size. `patchSize` can be a scalar (for
            %   square patches), or a two-element vector (for rectangular
            %   ones), in which case it's given in the order [nRows, nCols].
            %   If the image size is not a multiple integer of `patchSize`,
            %   a strip along the right and bottom edges of the image is
            %   ignored.
            %
            %   ImagePatchifier(image, patchSize, stride) uses the given
            %   stride for generating patches. Just like `patchSize`,
            %   `stride` can be a scalar or a two-element vector.
            %
            %   ImagePatchifier({imSize}, ...) generates image coordinates
            %   using only the image size, without reference to the actual
            %   image data. When using this format, calling the `get`
            %   function produces an error.
            
            if ~iscell(image)
                obj.image = image;
                obj.imSize = size(image);
            else
                obj.image = [];
                obj.imSize = image{1};
            end
            
            % handle scalar and two-element patch size
            if isscalar(patchSize)
                obj.patchSize = [patchSize patchSize];
            else
                obj.patchSize = patchSize;
            end
            
            % handle missing, scalar, and two-element stride
            if nargin < 3
                obj.stride = patchSize;
            elseif isscalar(stride)
                obj.stride = [stride stride];
            else
                obj.stride = stride;
            end
            
            obj.gridSize = 1 + floor((obj.imSize - obj.patchSize) ./ obj.stride);
            obj.gridIndex = [obj.gridSize(1) 0];
        end
        
        function ok = next(obj)
            % Go to the next location, return false if done.
                        
            if obj.gridIndex(1) < obj.gridSize(1)
                obj.gridIndex(1) = obj.gridIndex(1) + 1;
            else
                obj.gridIndex(1) = 1;
                obj.gridIndex(2) = obj.gridIndex(2) + 1;
            end
            ok = (obj.gridIndex(2) <= obj.gridSize(2));
        end
        
        function patch = get(obj)
            % Get the patch at the current location.

            if isempty(obj.image)
                error('ImagePatchifier:getnoimg', 'get requires actual image data.');
            end

            p = 1 + (obj.gridIndex - 1) .* obj.stride;
            patch = obj.image(p(1):p(1) + obj.patchSize(1) - 1, ...
                p(2):p(2) + obj.patchSize(2) - 1);
        end
        
        function n = numPatches(obj)
            % Return the number of patches.
            
            n = prod(obj.gridSize);
        end
        
        function coords = getPatchCoordinates(obj)
            % Get the coordinates for the current patch, in the form
            % [row1, col1, row2, col2].
            
            p = 1 + (obj.gridIndex - 1) .* obj.stride;
            coords = [p p+obj.patchSize-1];
        end
        
        function locs = patchLocations(obj)
            % Get matrix of patch locations.
            locStruct = generatePatchLocations(obj.imSize, obj.patchSize, obj.stride);
            locs = locStruct.patchLocations;
        end
        
        function idxs = gridIndices(obj)
            % Get matrix of grid indices.
            locStruct = generatePatchLocations(obj.imSize, obj.patchSize, obj.stride);
            idxs = locStruct.gridIndices;
        end
    end
    
end