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
    
    properties
        image = [];         % image we're working on
        patchSize = [];     % patch size we're using
        stride = [];        % stride we're using
%        numPatches = [];    % total number of patches
        gridSize = [];      % number of patches in each direction
        patchLocations = [];% numPatches x 2 matrix of patch locations in image
        gridIndices = [];   % numPatches x 2 matrix of grid indices
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
            
            obj.image = image;
            
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
            
            locStruct = generatePatchLocations(size(image), patchSize, obj.stride);
            obj.patchLocations = locStruct.patchLocations;
            obj.gridIndices = locStruct.gridIndices;
            obj.gridSize = locStruct.nPatches;
            obj.current_ = 0;
        end
        
        function ok = next(obj)
            % Go to the next location, return false if done.
                        
            if obj.current_ >= size(obj.patchLocations, 1)
                ok = false;
            else
                obj.current_ = obj.current_ + 1;
                obj.gridIndex = obj.gridIndices(obj.current_, :);
                ok = true;
            end
        end
        
        function patch = get(obj)
            % Get the patch at the current location.
            p = obj.patchLocations(obj.current_, :);
            patch = obj.image(p(1):p(1) + obj.patchSize(1) - 1, ...
                p(2):p(2) + obj.patchSize(2) - 1);
        end
        
        function n = numPatches(obj)
            % Return the number of patches.
            
            n = size(obj.patchLocations, 1);
        end
    end
    
end