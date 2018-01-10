classdef PatchAxisGenerator < handle
    % PatchAxisGenerator Generate texture patches along an axis in texture
    % space.
    %   
    %   obj = PatchAxisGenerator(group, axis, patchSize) creates a
    %   PatchAxisGenerator for a specific texture group, texture
    %   axis, and patch size. This allows to generate samples at a number
    %   of locations along that axis, ranging from completely unbiased
    %   patches at the center of the coordinate system, to the edge of the
    %   space.
    %
    %   Use the 'samples' method to generate samples at the current
    %   locations, and the 'next' method to move to the next location. Note
    %   that 'next' must be called once after object creation and before
    %   generating any samples. This initializes the generator.
    %
    %   The number of locations is set to a default value, but can be
    %   changed (using the 'nLocations' property) before the first call to
    %   'next'. Alternatively, the actual locations that are used can be
    %   set using the 'locations' property. A value of 0 corresponds to the
    %   center of the coordinate system, [1/n, ..., 1/n] (where 'n' is the
    %   number of gray levels), while 1 corresponds to the coordinates
    %   provided in the 'axis' argument to the constructor. Intermediate
    %   values are mapped to texture coordinates using linear interpolation.
    %
    %   The number of gray levels is inferred from the number of components
    %   in the texture axis. The texture axis is a probability distribution,
    %   so all the entries should be non-negative. It is automatically
    %   normalized so that the sum is 1 before texture generation.
    %
    %   The patch size can be a scalar or a 2-component vector, for
    %   generating rectangular patches.
    %
    %   Example:
    %
    %   % create generator for 3-gray-levels, 32x32 patches in the
    %   % [0.5 0.5 0] direction of the AB_1_2 texture plane
    %   generator = PatchAxisGenerator('AB_1_2', [0.5 0.5 0], 32);
    %
    %   % set the number of locations along the axis
    %   generator.nLocations = 16;
    %
    %   % choose number of samples per location
    %   nSamples = 8;
    %
    %   % generate
    %   while generator.next
    %       patches = generator.samples(nSamples);
    %   end
    
    properties
        nGray = [];         % number of gray levels
        patchSize = [];     % size of patch (scalar or 2-element vector)
        group = '';         % texture group
        axis = [];          % axis in texture space
        locations = [];     % locations where to sample
        nLocations = 16;    % number of locations to use
    end
    
    properties (Access=private)
        current_ = 0;
        btcDict_ = [];
        btcMakemapsOpts_ = [];
        auxOpts_ = [];
        mtcs_ = [];
        locations_ = [];
        tol_ = 1e-6;
    end
    
    methods
        function obj = PatchAxisGenerator(group0, axis0, patchSize0)
            if isscalar(patchSize0)
                patchSize0 = [patchSize0 patchSize0];
            end
            
            axis0 = axis0(:);
            obj.nGray = length(axis0);
            obj.patchSize = patchSize0;
            obj.group = group0;
            obj.axis = axis0/sum(axis0);
        end
        
        function init(obj)
            % Initialize the texture generator.
            
            % overall map setup
            obj.btcDict_ = btc_define;
            %aug_opts = [];
            %aug_opts.ifstd = 1;
            obj.btcMakemapsOpts_ = [];
            obj.btcMakemapsOpts_.nmaps = 1;
            obj.btcMakemapsOpts_.show = 0;
            obj.btcMakemapsOpts_.area = obj.patchSize;
            obj.auxOpts_ = btc_auxopts; % won't need metropolis arguments here
            
            % 3 grayscale levels
            obj.mtcs_ = mtc_define(obj.nGray);
            
            % set up the locations, unless explicitly provided
            if ~isempty(obj.locations)
                obj.locations_ = obj.locations;
            else
                % 1/n*(1-t) + a*t >= 0
                % 1/n + (a-1/n)*t >= 0
                %   a) if a>1/n, then t >= -1/n / (a-1/n) = -1/(a*n - 1)
                %       --> not a problem
                %   b) if a<1/n, then t <= -1/n / (a-1/n) = 1/(1 - a*n)
                %       --> t <= min(1/(1-a*n)) over all axes
                %           t <= 1/max(1-a*n)
                %           t <= 1/(1 - min(a)*n)
                minAxisCoord = min(obj.axis);
                if abs(minAxisCoord*obj.nGray - 1) > obj.tol_
                    maxT = 1/(1 - minAxisCoord*obj.nGray);
                    obj.locations_ = linspace(0, maxT, obj.nLocations);
                else
                    obj.locations_ = zeros(1, obj.nLocations);
                end
            end
            % start at the first location
            obj.current_ = 0;
        end
        
        function ok = next(obj)
            % Go to the next location, return false if done.
            
            % initialize object if needed
            if obj.current_ == 0
                obj.init
            end
            
            if obj.current_ >= length(obj.locations_)
                ok = false;
            else
                obj.current_ = obj.current_ + 1;
                ok = true;
            end
        end
        
        function res = samples(obj, n)
            % Generate texture samples.
            %   obj.SAMPLES generates ones patch.
            %   obj.SAMPLES(n) generates 'n' patches, arranged in a 3d
            %   array.
            if nargin < 2
                n = 1;
            end
            
            if obj.current_ < 1
                error('PatchAxisGenerator:sample', 'The ''next'' method must be called before samples can be generated.');
            end
            
            res = zeros([obj.patchSize(:)' n]);
            crtLocation = obj.locations_(obj.current_);
            crtCoords = ones(1, obj.nGray)/obj.nGray*(1 - crtLocation) + ...
                obj.axis(:)'*crtLocation;
            
            % set up the texture generator for this particular direction
            cgsStruct = [];
            cgsStruct.(obj.group) = crtCoords;
            augCoords = mtc_augcoords(cgsStruct, obj.mtcs_, obj.btcDict_);

            for i = 1:n
                res(:, :, i) = btc_makemaps(augCoords.method{1}, ...
                    obj.btcMakemapsOpts_, obj.btcDict_, obj.auxOpts_, []) / ...
                    (obj.nGray - 1);
            end
        end
    end
    
end