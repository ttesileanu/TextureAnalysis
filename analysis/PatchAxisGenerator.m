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
    %   More specifically, the texture axis is the line connecting the
    %   origin of texture space, [1/n, ..., 1/n] (where 'n' is the number
    %   of gray levels), and the location identified by the 'axis'
    %   argument.
    %
    %   Use the 'samples' method to generate samples at the current
    %   locations, and the 'next' method to move to the next location. Note
    %   that 'next' must be called once after object creation and before
    %   generating any samples. This initializes the generator.
    %
    %   The number of locations is set to a default value, but can be
    %   changed (using the 'nLocations' property) before the first call to
    %   'next'. The start and end locations are set to 0 and the maximum
    %   allowed value, respectively, but can be changed using the
    %   'minLocation' and 'maxLocation' properties. Alternatively, the
    %   actual locations that are used can be set using the 'locations'
    %   property. A value of 0 corresponds to the center of the coordinate
    %   system, 1/n, ..., 1/n], while 1 corresponds to the coordinates
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
    %   obj = PatchAxisGenerator({group1, group2, ...}, axesMatrix, patchSize)
    %   generates patches along an axis in which several different texture
    %   groups are set. The directions along which the samples are
    %   generated in each texture group are given by the rows of the axis
    %   parameter.
    %
    %   Example single group:
    %
    %       % create generator for 3 gray-levels, 32x32 patches in the
    %       % [0.5 0.5 0] direction of the AB_1_2 texture plane
    %       generator = PatchAxisGenerator('AB_1_2', [0.5 0.5 0], 32);
    %
    %       % set the number of locations along the axis
    %       generator.nLocations = 16;
    %
    %       % choose number of samples per location
    %       nSamples = 8;
    %
    %       % generate
    %       while generator.next
    %           patches = generator.samples(nSamples);
    %       end
    %
    %   Example multiple groups:
    %
    %       % create generator for 2 gray-levels, 32x64 patches in the
    %       % [0.9 0.1] direction of the AD_1_1 texture plane and the
    %       % [0 1] direction of the ABC_1_1_1 texture plane
    %       generator = PatchAxisGenerator({'AD_1_1', 'ABC_1_1_1'}, ...
    %           [0.9 0.1 ; 0 1], [32 64]);
    %
    %       % set the number of locations along the axis
    %       generator.nLocations = 10;
    %
    %       % generate one patch at a time
    %       while generator.next
    %           patch = generator.samples;
    %       end
    
    properties
        nGray = [];         % number of gray levels
        patchSize = [];     % size of patch (scalar or 2-element vector)
        group = '';         % texture group(s)
        axis = [];          % axis/axes in texture space
        locations = [];     % locations where to sample
        nLocations = 16;    % number of locations to use
        minLocation = 0;    % starting location
        maxLocation = [];   % ending location
        mcIterations = 100; % number of iterations to use for Monte Carlo
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
            % handle square or rectangular patch sizes
            if isscalar(patchSize0)
                patchSize0 = [patchSize0 patchSize0];
            end
            % handle single or multi-group
            if ~iscell(group0)
                group0 = {group0};
            end
            % make sure we have an axis for each texture group
            if length(group0) ~= size(axis0, 1)
                error('PatchAxisGenerator:badngrp', ...
                    'The number of rows of the axis argument should match the number of groups.');
            end
            
            % store the parameters
            obj.nGray = size(axis0, 2);
            obj.patchSize = patchSize0;
            obj.group = group0;
            
            % make sure that each axis is a properly normalized probability
            % distribution
%             if any(axis0(:) < 0)
%                 error('PatchAxisGenerator:negax', 'All axis components should be non-negative.');
%             end
            obj.axis = bsxfun(@rdivide, axis0, sum(axis0, 2));
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
            obj.auxOpts_ = btc_auxopts;
            obj.auxOpts_.metro_opts.numiters = obj.mcIterations;
            
            % set number of grayscale levels
            obj.mtcs_ = mtc_define(obj.nGray);
            
            % set up the locations, unless explicitly provided
            if ~isempty(obj.locations)
                obj.locations_ = obj.locations;
            else
                if isempty(obj.maxLocation)
                    if length(obj.group) == 1
                        % 1/n*(1-t) + a*t >= 0
                        % 1/n + (a-1/n)*t >= 0
                        %   a) if a>1/n, then t >= -1/n / (a-1/n) = -1/(a*n - 1)
                        %       --> not a problem; can't happen for all coords
                        %   b) if a<1/n, then t <= -1/n / (a-1/n) = 1/(1 - a*n)
                        %       --> t <= min(1/(1-a*n)) over all axes
                        %           t <= 1/max(1-a*n)
                        %           t <= 1/(1 - min(a)*n)
                        %   c) if all coordinates = 1/n, then we're just staying at
                        %       the center
                        minAxisCoord = min(obj.axis(:));
                        if abs(minAxisCoord*obj.nGray - 1) > obj.tol_
                            maxT = 1/(1 - minAxisCoord*obj.nGray);
                        else
                            maxT = 0;
                        end
                    else
                        % bracketing between 0 and 1: see what works
                        obj.current_ = 1;
                        % first check the simple cases
                        obj.locations_ = 1;
                        if ~isempty(obj.samples('check'))
                            maxT = 1;
                        else
                            l = 0;
                            r = 1;
                            % invariants: l always works, r always doesn't
                            while r - l > obj.tol_
                                m = (l + r) / 2;
                                obj.locations_ = m;
                                if isempty(obj.samples('check'))
                                    r = m;
                                else
                                    l = m;
                                end
                            end
                            % we pinpointed point where generation doesn't
                            % work to within obj.tol_
                            maxT = l;
                        end
                    end
                else
                    maxT = obj.maxLocation;
                end
                obj.locations_ = linspace(0, maxT, obj.nLocations);
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
                error('PatchAxisGenerator:noinit', ...
                    'The ''next'' method must be called once before samples can be generated.');
            end
            
            crtLocation = obj.locations_(obj.current_);
            crtCoords = ones(size(obj.axis))/obj.nGray*(1 - crtLocation) + ...
                obj.axis*crtLocation;
            
            % set up the texture generator for this particular direction
            cgsStruct = [];
            for k = 1:length(obj.group)
                cgsStruct.(obj.group{k}) = crtCoords(k, :);
            end
            augCoords = mtc_augcoords(cgsStruct, obj.mtcs_, obj.btcDict_);
            
            if isempty(augCoords.method)
                res = [];
                return;
            end
            
            if strcmp(n, 'check')
                % just check that we could generate samples in principle
                res = 1;
                return;
            end

            res = zeros([obj.patchSize(:)' n]);
            for i = 1:n
                res(:, :, i) = btc_makemaps(augCoords.method{1}, ...
                    obj.btcMakemapsOpts_, obj.btcDict_, obj.auxOpts_, []) / ...
                    (obj.nGray - 1);
            end
        end
        
        function locs = getLocations(obj)
            locs = obj.locations_;
        end
    end
    
end