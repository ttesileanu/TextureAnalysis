classdef MatrixPlotter < handle
    % MatrixPlotter Arrange a series of plots in a grid.
    %
    %   This class simplifies the creation of a figure containing many
    %   plots by arranging them in a grid. It automatically chooses an
    %   appropriate size for the figure and selects a balanced number of
    %   rows and columne.
    %
    %   Example usage:
    %
    %       plotter = MatrixPlotter(13); % 13 plots
    %       while plotter.next
    %           scatter(randn(100, 1), randn(100, 1), '.k');
    %           title(int2str(plotter.index));
    %       end
    
    properties
        fig_aspect = 3/2;   % approximate aspect ratio of figure
        ax_aspect = 3/2;    % approximate aspect ratio of each axes object
        screen_edge = 0.15;  % fraction of sreen that should be left unused
    end
    
    properties (SetAccess=private)
        index = [];    
    end
    
    properties (SetAccess=immutable)
        n_plots = [];
    end
    
    properties (Access=private)
        fig_ = [];
        ax_ = [];
        nx_ = [];
        ny_ = [];
    end
    
    methods
        function obj = MatrixPlotter(n_plots)
            % MatrixPlotter Constructor.
            %   MatrixPlotter(n_plots) creates the figure, sized for the 
            %   given number of plots, and initializes the MatrixPlotter
            %   object.
            
            obj.n_plots = n_plots;
            obj.createFigure_;
        end
        
        function good = next(obj)
            % NEXT Create next set of axes.
            %   The function returns `true` if the access were created,
            %   `false` if all the axes have already been generated.
            
            % check if we've already created all the plots
            if obj.index >= obj.n_plots
                good = false;
                return;
            end
            
            % generate the next axes
            ax = axes;
            px = mod(obj.index, obj.nx_);
            py = floor(obj.index/obj.nx_);
            
            ax.OuterPosition = [px/obj.nx_ 1-(py+1)/obj.ny_ 1/obj.nx_ 1/obj.ny_];

            % store the axes and return
            obj.index = obj.index + 1;
            obj.ax_(obj.index) = ax;
            
            good = true;
        end
        
        function ax = get_axes(obj)
            % Get a copy of the axes objects.
            ax = obj.ax_;
        end
        
        function fig = get_figure(obj)
            % Get a copy of the figure object.
            fig = obj.fig_;
        end
    end
    
    methods (Access=private)
        function createFigure_(obj)
            % Calculate appropriate size and create the figure.
            
            % figure out how to place plots
            obj.nx_ = ceil(sqrt(obj.fig_aspect/obj.ax_aspect*obj.n_plots));
            obj.ny_ = ceil(obj.n_plots / obj.nx_);

            % figure out maximum figure size
            screenRects = get(groot, 'ScreenSize');
            screen_x = screenRects(3);
            screen_y = screenRects(4);
            max_x = screen_x*(1 - obj.screen_edge);
            max_y = screen_y*(1 - obj.screen_edge);
            
            % figure out the size of the figure and of each axis
            fig = figure;
            ax_x = max_x/obj.nx_;
            ax_y = ax_x/obj.ax_aspect;
            if ax_y > max_y/obj.ny_
                ax_y = max_y/obj.ny_;
                ax_x = ax_y*obj.ax_aspect;
            end
            fig_x = ax_x*obj.nx_;
            fig_y = ax_y*obj.ny_;
            edge_x = (screen_x - fig_x)/2;
            edge_y = (screen_y - fig_y)/2;
            fig.Position = [edge_x edge_y fig_x fig_y];
            
            obj.fig_ = fig;
            
            % initialize the axis iteration
            obj.index = 0;
        end
    end
    
end