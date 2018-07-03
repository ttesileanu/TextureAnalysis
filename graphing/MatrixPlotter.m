classdef MatrixPlotter < handle
    % MatrixPlotter Arrange a series of plots in a grid.
    %
    %   This class simplifies the creation of a figure containing many
    %   plots by arranging them in a grid. By default, it automatically
    %   chooses an appropriate size for the figure and selects a balanced
    %   number of rows and columns.
    %
    %   With additional options, the class can make figures with more
    %   complex layouts. See the example usage below.
    %
    %   Examples of usage:
    %
    %   1. Basic usage
    %
    %       plotter = MatrixPlotter(13); % 13 plots
    %       while plotter.next
    %           scatter(randn(100, 1), randn(100, 1), '.k');
    %           title(int2str(plotter.index));
    %       end
    %
    %   2. Complex layout
    %
    %       plotter = MatrixPlotter(10, 'fixedShape', [5 2]);
    %       plotter.fixSize([10 4], 'inches');
    %       while plotter.next
    %           hist(randn(100, 1));
    %       end
    %
    %   3. Using external axes
    %
    %       fig = figure;
    %       ax = zeros(1, 2);
    %       ax(1) = axes('OuterPosition', [0 0 0.5 0.5]);
    %       ax(2) = axes('OuterPosition', [0.5 0.5 0.5 0.5]);
    %       plotter = MatrixPlotter(2, 'fixedAxes', ax);
    %       while plotter.next
    %           plot(randn(100, 1), randn(100, 1));
    %       end
    
    properties
        figAspect = 3/2;    % approximate aspect ratio of figure
        axAspect = 3/2;     % approximate aspect ratio of each axes object
        screenEdge = 0.15;  % fraction of sreen that should be left unused
        fixedShape = [];    % force shape, [nrows, ncols]
        fixedSize = [];     % force a size for the figure, [figx, figy]
        fixedSizeUnits = [];% units to use to force size
        fixedAxes = [];     % force the use of a given set of axes
    end
    
    properties (SetAccess=private)
        index = [];    
    end
    
    properties (SetAccess=immutable)
        nPlots = [];
    end
    
    properties (Access=private)
        fig_ = [];
        ax_ = [];
        nx_ = [];
        ny_ = [];
    end
    
    methods
        function obj = MatrixPlotter(nPlots, varargin)
            % MatrixPlotter Constructor.
            %   MatrixPlotter(nPlots) creates the figure, sized for the 
            %   given number of plots, and initializes the MatrixPlotter
            %   object.
            
            parser = inputParser;
            parser.CaseSensitive = true;
            parser.FunctionName = [mfilename '.constructor'];
            
            none = {'none'};
            parser.addParameter('figAspect', none, @(x) isnumeric(x) && isscalar(x) && x > 0);
            parser.addParameter('axAspect', none, @(x) isnumeric(x) && isscalar(x) && x > 0);
            parser.addParameter('screenEdge', none, @(x) isnumeric(x) && isscalar(x) && x >= 0);
            parser.addParameter('fixedShape', none, @(v) isnumeric(v) && isvector(v) && length(v) == 2 && all(v > 0));
            parser.addParameter('fixedSize', none, @(v) isnumeric(v) && isvector(v) && length(v) == 2 && all(v > 0));
            parser.addParameter('fixedSizeUnits', none, @(s) ischar(s) && isvector(s));
            parser.addParameter('fixedAxes', none, @(v) isvector(v));
            
            % parse
            parser.parse(varargin{:});
            params = parser.Results;
            
            % update the properties
            fields = fieldnames(params);
            for i = 1:numel(fields)
                crtValue = params.(fields{i});
                if ~isequal(crtValue, none)
                    obj.(fields{i}) = crtValue;
                end
            end
            
            obj.nPlots = nPlots;
        end
        
        function fixShape(obj, shape)
            % obj.fixShape(shape)
            %   Set the shape of the figure in number of plots ([rows, cols]).
            obj.fixedShape = shape;
        end
        
        function fixSize(obj, size, sizeUnits)
            % obj.fixShape(size)
            %   Set the display size of the figure ([x, y] order!).
            %
            % obj.fixShape(size, sizeUnits)
            %   Set the size in the given units.
            
            if nargin >= 2
                obj.fixedSizeUnits = sizeUnits;
            end
            obj.fixedSize = size;
        end
        
        function fixAxes(obj, ax)
            % obj.fixAxes(ax)
            %   Use the given array of axes for the plots. In this case no
            %   new figure is created and the pre-existing axes are used.
            
            obj.fixedAxes = ax;
        end
        
        function good = next(obj)
            % NEXT Create next set of axes.
            %   The function returns `true` if the access were created,
            %   `false` if all the axes have already been generated.
            
            % make sure we have a figure
            if isempty(obj.fig_)
                obj.createFigure_;
            end
            
            % check if we've already created all the plots
            if obj.index >= obj.nPlots
                good = false;
                return;
            end
            
            if isempty(obj.fixedAxes)
                % generate the next axes
                ax = axes;
                px = mod(obj.index, obj.nx_);
                py = floor(obj.index/obj.nx_);
                
                ax.OuterPosition = [px/obj.nx_ 1-(py+1)/obj.ny_ 1/obj.nx_ 1/obj.ny_];
                
                % store the axes and return
                obj.index = obj.index + 1;
                obj.ax_(obj.index) = ax;
            else
                % activate the proper axis
                obj.index = obj.index + 1;
                axes(obj.ax_(obj.index));
            end
            
            good = true;
        end
        
        function ax = getAxes(obj)
            % Get a copy of the axes objects.
            ax = obj.ax_;
        end
        
        function fig = getFigure(obj)
            % Get a copy of the figure object.
            fig = obj.fig_;
        end
    end
    
    methods (Access=private)
        function createFigure_(obj)
            % Calculate appropriate size and create the figure.
            
            % we generate our own figure
            if isempty(obj.fixedAxes)
                if isempty(obj.fixedShape)
                    % figure out how to place plots
                    obj.nx_ = ceil(sqrt(obj.figAspect/obj.axAspect*obj.nPlots));
                    obj.ny_ = ceil(obj.nPlots / obj.nx_);
                else
                    obj.nx_ = obj.fixedShape(2);
                    obj.ny_ = obj.fixedShape(1);
                end
                
                % figure out maximum figure size
                screenRects = get(groot, 'ScreenSize');
                screenX = screenRects(3);
                screenY = screenRects(4);
                
                fig = figure;
                if isempty(obj.fixedSize)
                    maxX = screenX*(1 - obj.screenEdge);
                    maxY = screenY*(1 - obj.screenEdge);
                    
                    % figure out the size of the figure and of each axis
                    axX = maxX/obj.nx_;
                    axY = axX/obj.axAspect;
                    if axY > maxY/obj.ny_
                        axY = maxY/obj.ny_;
                        axX = axY*obj.axAspect;
                    end
                    figX = axX*obj.nx_;
                    figY = axY*obj.ny_;
                else
                    figX = obj.fixedSize(1);
                    figY = obj.fixedSize(2);
                    
                    if ~isempty(obj.fixedSizeUnits)
                        old_units = fig.Units;
                        fig.Units = obj.fixedSizeUnits;
                        fig.Position = [0 0 figX figY];
                        fig.Units = old_units;
                        figX = fig.Position(3);
                        figY = fig.Position(4);
                    end
                end
                edgeX = (screenX - figX)/2;
                edgeY = (screenY - figY)/2;
                fig.Position = [edgeX edgeY figX figY];
                
                obj.fig_ = fig;
            else
                obj.fig_ = gcf;
                obj.ax_ = obj.fixedAxes;
            end
            
            % initialize the axis iteration
            obj.index = 0;
        end
    end
    
end