classdef TextProgress < handle
    % TextProgress Text-mode progress bar.
    %
    %   Based on textprogressbar entry on MATLAB Central:
    %   https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar
    %
    %   Usage examples:
    %
    %   1. Progress bar without any other text.
    %
    %       progress1 = TextProgress;
    %       for i = 1:10
    %           progress1.update(10*i);
    %           pause(0.5);
    %       end
    %       progress1.done;     % needed for newline
    %
    %   2. Progress bar with prefix and suffix, custom length
    %
    %       progress2 = TextProgress('working...');
    %       progress2.length = 25;
    %       for i = 1:100
    %           progress2.update(i);
    %           pause(0.05);
    %       end
    %       progress2.done('done!');
    %
    %   3. Progress bar with changing prefix, suffix, fixed prefix spacing
    %
    %       progress3 = TextProgress('', 'prespace', 24);
    %       words = {'This', 'is', 'a', 'test'};
    %       for i = 1:length(words)
    %           progress3.update(i*100/length(words), 'prefix', words{i}, ...
    %               'suffix', ['word #' int2str(i)]);
    %           pause(1);
    %       end
    %       progress3.done('done!', 'prefix', 'done!');
    %
    %   4. Change the update interval
    %
    %       progress4 = TextProgress('test', 'updinterval', 1.0);
    %       for i = 1:100
    %           progress4.update(i);
    %           pause(0.1);
    %       end
    %       progress4.done;     % needed for newline

    
    properties
        length = 10;
        prespace = [];
        prefix = '';
        suffix = '';
        percent = 0;
        timer = true;
        timefmt = '%5.1f';
        updinterval = 0.25;
    end
    
    properties (Access=private)
        backStr_ = '';
        t0_ = [];
        tlastupd_ = [];
    end
    
    methods
        function obj = TextProgress(varargin)
            % TextProgress Constructor.
            %   TextProgress initializes the progress bar with default
            %   options and no prefix.
            %
            %   TextProgress('text') uses the given text as prefix.
            %
            %   TextProgress(..., 'key', value, ...) sets the options after
            %   initializing the bar.
            if isempty(varargin)
                obj.update_;
            else
                obj.update_('prefix', varargin{:});
            end
            obj.t0_ = tic;
        end
        
        function update(obj, percent, varargin)
            % UPDATE Update the progress bar.
            %   UPDATE(percent) updates the percentage value for the bar
            %   (`percent` must be a numeric value).
            %
            %   UPDATE('key', value, ...) updates the bar by changing the
            %   options.
            %
            %   Additionally, the option 'force', true can be used to force
            %   update even if the appropriate interval since the last
            %   update hasn't yet passed.
            if isnumeric(percent)
                obj.update_('percent', percent, varargin{:});
            else
                obj.update_(percent, varargin{:});
            end
        end
        
        function done(obj, suffix, varargin)
            % DONE Finish with the progress bar.
            %   DONE Sets the percentage to 100 and displays a newline.
            %
            %   DONE('text') also sets the suffix to the given text.
            %
            %   DONE(..., 'key', value, ...) also sets other options.
            if nargin >= 2
                obj.update_('percent', 100, 'suffix', suffix, 'force', true, ...
                    varargin{:});
            else
                obj.update_('percent', 100, 'force', true);
            end
            fprintf('\n');
            obj.backStr_ = '';
        end
        
        function update_(obj, varargin)
            % Update the progress bar using the given key/value pairs.
            % parse optional arguments
            parser = inputParser;
            parser.CaseSensitive = true;
            parser.FunctionName = [mfilename '.update_'];
            
            none = {'none'};
            parser.addParameter('length', none, @(n) isnumeric(n) && isreal(n) && n > 3 && floor(n) == n);
            parser.addParameter('prespace', none, @(n) isnumeric(n) && isreal(n) && n > 3 && floor(n) == n);
            parser.addParameter('prefix', none, @(s) ischar(s) && (isvector(s) || isempty(s)));
            parser.addParameter('suffix', none, @(s) ischar(s) && (isvector(s) || isempty(s)));
            parser.addParameter('percent', none, @(n) isnumeric(n) && isreal(n));
            parser.addParameter('timer', none, @(b) islogical(b) && isscalar(b));
            parser.addParameter('timefmt', none, @(s) ischar(s) && (isvector(s) || isempty(s)));
            parser.addParameter('updinterval', none, @(n) isnumeric(n) && isreal(n));
            parser.addParameter('force', false, @(b) islogical(b) && isscalar(b));
            
            % parse
            parser.parse(varargin{:});
            params = parser.Results;

            % force is only passed to display_
            force = params.force;
            params = rmfield(params, 'force');
            
            % update the properties
            fields = fieldnames(params);
            for i = 1:numel(fields)
                crt_value = params.(fields{i});
                if ~isequal(crt_value, none)
                    obj.(fields{i}) = crt_value;
                end
            end

            % handle the display
            obj.display_(force);
        end
        
        function display_(obj, force)
            % Display the bar. Normally this is only done if enough time
            % passed since the last update. This can be overridden by
            % setting `force` to `true`.
            if nargin < 2
                force = false;
            end
            
            if ~force && ~isempty(obj.tlastupd_) && toc(obj.tlastupd_) < obj.updinterval
                % don't update too often
                return;
            end
            if isempty(obj.prespace)
                bar = obj.prefix;
            else
                bar = sprintf(['%-' int2str(obj.prespace) 's'], obj.prefix);
            end
            
            if ~isempty(bar)
                bar = [bar ' '];
            end
            obj.percent = min(max(obj.percent, 0), 100);
            bar = [bar sprintf('%3d%%', round(obj.percent)) ' ['];
            nDots = round(obj.percent*obj.length/100);
            bar = [bar repmat('.', 1, nDots) repmat(' ', 1, obj.length - nDots)];
            bar = [bar ']'];
            
            if obj.timer
                dt = 0;
                if ~isempty(obj.t0_)
                    dt = toc(obj.t0_);
                end
                bar = [bar ' ' sprintf(obj.timefmt, dt) ' sec'];
            end
            
            if ~isempty(obj.suffix)
                bar = [bar ' ' obj.suffix];
            end
            % add an extra space to separate from prompt
            bar = [bar ' '];
            
            fprintf(1, [obj.backStr_ '%s'], bar);
            obj.backStr_ = repmat('\b', 1, numel(bar));
            
            obj.tlastupd_ = tic;
        end
    end
    
end
