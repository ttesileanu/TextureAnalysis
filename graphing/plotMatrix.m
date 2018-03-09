function plotMatrix(m, varargin)
% plotMatrix Plot the entries of a matrix.
%   plotMatrix(m) makes a figure showing the contents of the matrix `m`.
%   This behaves like Matlab's imagesc, except for small matrices it draws
%   a grid of lines in-between the entries. It can also display the values
%   for each entry.
%
%   plotMatrix(m, [clow chigh]) scales the entries so that clow is mapped
%   to the first entry in the colormap, and chigh is mapped to the last
%   entry in the colormap.
%
%   Options:
%    'edge': number
%       Amount of space to leave between matrix entries, as a fraction of
%       the space available for each entry.
%    'frameColor': Matlab color specification, or empty
%       If not empty, this should be a letter or an RGB triplet giving the
%       color of the frame surrounding the plot. If empty, no frame is
%       drawn. By default, there's a white frame whenever the image has 32
%       or fewer rows and 32 or fewer columns.
%    'fontSize': number
%       Size of font used to display values of entries as a fraction of
%       each square in the drawing. Set to an empty matrix to not display
%       entries (default), or 'auto' to use the default size.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('clim', [], @(v) isvector(v) && isnumeric(v) && length(v) == 2);

parser.addParameter('edge', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('frameColor', 'auto');
parser.addParameter('fontSize', [], @(x) isempty(x) || strcmp(x, 'auto') || ...
    isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check inputs
if ~ismatrix(m) || ~isnumeric(m)
    error([mfilename ':badshape'], 'Input should be a numeric matrix.');
end

% handle some defaults
if strcmp(params.fontSize, 'auto')
    params.fontSize = 0.3;
end

if strcmp(params.frameColor, 'auto')
    if all(size(m) <= 32)
        params.frameColor = [1 1 1];
    else
        params.frameColor = [];
    end
end

% draw matrix entries
imagesc(m, params.clim);
hold on;

% add grid
if ~isempty(params.frameColor)
    [M, N] = size(m);
    
    % convert units for params.edge
    ax = gca;
    oldUnits = ax.Units;
    ax.Units = 'points';
    imPos = ax.Position;
    ax.Units = oldUnits;
    
    width = params.edge * imPos(3) / N;
    for i = 0:M
        line([0.5, N + 0.5], i + [0.5, 0.5], 'color', params.frameColor, 'linewidth', width);
    end
    for i = 0:N
        line(i + [0.5, 0.5], [0.5, M + 0.5], 'color', params.frameColor, 'linewidth', width);
    end
end

% add text, if required
% XXX implement

% make axis square, and remove box
axis image;
%box off;
%axis off;

% make fonts nicer
beautifygraph;

end