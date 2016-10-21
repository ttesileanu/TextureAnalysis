function drawSymmetricMatrix(m, varargin)
% drawSymmetricMatrix Make a diagram of one half of a symmetric matrix.
%   drawSymmetricMatrix(m) makes a figure showing the contents of the
%   symmetric matrix `m`. Only one half of the matrix is shown (the sub-
%   diagonal one by default). The current colormap is used, and the colors
%   will be updated automatically if the colormap is changed afterwards.
%
%   drawSymmetricMatrix(m, [clow chigh]) scales the entries so that clow is
%   mapped to the first entry in the colormap, and chigh is mapped to the
%   last entry in the colormap.
%
%   Options:
%    'abs': bool
%       If true, choose color based on absolute value of entries of `m`
%       (but text still uses the actual value).
%    'edge': number
%       Amount of space to leave between matrix entries, as a fraction of
%       the space available for each entry.
%    'frameColor': Matlab color specification, or empty
%       If not empty, this should be a letter or an RGB triplet giving the
%       color of the frame surrounding the plot. If empty, no frame is
%       drawn.
%    'half': string
%       Identify the half of the matrix to show, 'upper' or 'lower'.
%    'fontSize': number
%       Size of font used to display values of entries as a fraction of
%       each square in the drawing. Set to an empty matrix to not display
%       entries (default), or 'auto' to use the default size.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('clim', [], @(v) isvector(v) && isnumeric(v) && length(v) == 2);

parser.addParameter('abs', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('edge', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('frameColor', [0 0 0]);
parser.addParameter('half', 'lower', @(s) ismember(s, {'lower', 'upper'}));
parser.addParameter('fontSize', [], @(x) isempty(x) || strcmp(x, 'auto') || ...
    isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check inputs
n = size(m, 1);
if ~ismatrix(m) || ~isnumeric(m) || size(m, 2) ~= n
    error([mfilename ':badshape'], 'Input should be a square numeric matrix.');
end

% handle some defaults
if strcmp(params.fontSize, 'auto')
    params.fontSize = 0.3;
end

%colormap(gca, 'gray');

% draw the heatmap
% fill with NaNs the side that we don't want drawn
if strcmp(params.half, 'lower')
    m = m + triu(nan(n), 1);
else
    m = m + tril(nan(n), -1);
end
hEdge = params.edge/2;
for i = 1:n
    for j = 1:n
        crt = m(i, j);
        if isnan(crt)
            continue;
        end
        if params.abs
            crtCol = abs(crt);
        else
            crtCol = crt;
        end
        patch([j-0.5+hEdge j+0.5-hEdge j+0.5-hEdge j-0.5+hEdge], ...
              [i-0.5+hEdge i-0.5+hEdge i+0.5-hEdge i+0.5-hEdge], crtCol, ...
            'edgecolor', 'none');
        
        % write in the numbers
        if ~isempty(params.fontSize)
            txt = num2str(crt, '%.2f');
            text(j, i, txt, 'horizontalalignment', 'center', ...
                'verticalalignment', 'middle', ...
                'fontunits', 'normalized', 'fontsize', params.fontSize/n);
        end
    end
end

% get rid of axes
axis('image');
axis('off');

% reverse y axis
set(gca, 'ydir', 'reverse');

% set color range
if ~isempty(params.clim)
    set(gca, 'clim', params.clim);
end

% set view range
xlim([0.5 n+0.5]);
ylim([0.5 n+0.5]);

% draw the frame
if ~isempty(params.frameColor)
    hold on;
    
    minP = 0.5 - hEdge;
    maxP = n + 0.5 + hEdge;
    
    x1 = [minP-params.edge minP minP maxP+params.edge maxP+params.edge minP-params.edge];
    y1 = [minP-params.edge minP-params.edge maxP maxP maxP+params.edge maxP+params.edge];
    if strcmp(params.half, 'upper')
        x1 = n + 1 - x1;
        y1 = n + 1 - y1;
    end
    fill(x1, y1, params.frameColor, 'edgecolor', 'none');
    
    pts2 = zeros(4*n+2, 2);
    pts2(1, :) = [minP-params.edge minP-params.edge];
    for i = 1:n
        pts2(2*i, :) = [minP+2*params.edge+i minP-params.edge+(i-1)];
        pts2(2*i+1, :) = [minP+2*params.edge+i minP-params.edge+i];
    end
    % lengthening last segment a bit
    pts2(2*n+1, :) = [maxP+params.edge maxP+params.edge];
    pts2(2*n+2, :) = [maxP maxP+params.edge];
    for i = 1:n
        pts2(2*n+2*i+1, :) = [maxP-i+1 maxP-i-params.edge];
        pts2(2*n+2*i+2, :) = [maxP-i maxP-i-params.edge];
    end
    pts2(4*n+2, :) = [minP-params.edge minP];
    if strcmp(params.half, 'upper')
        pts2 = n + 1 - pts2;
    end
    fill(pts2(:, 1), pts2(:, 2), params.frameColor, 'edgecolor', 'none');
    
    xlim([minP - params.edge, maxP+params.edge]);
    ylim([minP - params.edge, maxP+params.edge]);
end

end