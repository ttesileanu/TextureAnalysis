function drawStochasticMatrix(m, varargin)
% drawStochasticMatrix Make a diagram of a stochastic matrix.
%   drawStochasticMatrix(m) makes a figure showing the contents of the
%   stochastic matrix `m`. Each entry is colored according to a blue -
%   white - red color map where red = 1 and white = uniform distribution.
%   The probability values are also overlayed on top of the colored
%   squares.
%
%   Options:
%    'edge': number
%       Amount of space to leave between matrix entries, as a fraction of
%       the space available for each entry.
%    'edgeColor': Matlab color specification
%       This can be a letter or an RGB triplet giving the color of the
%       frame surrounding the matrix entries.
%    'fontSize': number
%       Size of font used to display values of entries as a fraction of
%       each square in the drawing.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('edge', 0.1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('edgeColor', [0 0 0]);
parser.addParameter('fontSize', 0.3, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check inputs
n = size(m, 1);
if ~ismatrix(m) || ~isnumeric(m) || size(m, 2) ~= n
    error([mfilename ':badshape'], 'Input should be a square numeric matrix.');
end

% set the colormap
nLevels = 255;
cmap = zeros(nLevels, 3);
midPoint = floor((nLevels+1)/2);
gradient = linspace(0, 1, midPoint-1)';
cmap(midPoint, :) = [1 1 1];
cmap(1:midPoint-1, :) = [gradient gradient ones(midPoint-1, 1)];
cmap(midPoint+1:end, :) = [ones(midPoint-1, 1) flipud(gradient) flipud(gradient)];

colormap(gca, cmap);

% draw the heatmap
% set 1/n in the middle (white), make 1 full red; 0 will be a shade of blue
% depending on n (for n = 2, it would be full blue)
imagesc(m, [1-2*(n-1)/n 1]);

% get rid of axes
axis('image');
axis('off');

% draw the frame
hEdge = params.edge/2;
minP = 0.5-hEdge;
maxP = n+0.5+hEdge;
xlim([minP maxP]);
ylim([minP maxP]);

hold on;
fill([minP minP+params.edge minP+params.edge minP], [minP minP maxP maxP], params.edgeColor, ...
    'edgecolor', 'none');
fill([maxP maxP-params.edge maxP-params.edge maxP], [minP minP maxP maxP], params.edgeColor, ...
    'edgecolor', 'none');
fill([minP minP maxP maxP], [maxP maxP-params.edge maxP-params.edge maxP], params.edgeColor, ...
    'edgecolor', 'none');
fill([minP minP maxP maxP], [minP minP+params.edge minP+params.edge minP], params.edgeColor, ...
    'edgecolor', 'none');
for i = 1:n-1
    fill([minP+i minP+i+params.edge minP+i+params.edge minP+i], ...
        [minP minP maxP maxP], params.edgeColor, 'edgecolor', 'none');
    fill([minP minP maxP maxP], ...
        [minP+i minP+i+params.edge minP+i+params.edge minP+i], ...
        params.edgeColor, 'edgecolor', 'none');
end

% write in the numbers
for i = 1:n
    for j = 1:n
        proc = max(min(round(m(i, j)*100), 100), 0);
        if proc == 100
            txt = '1.0';
        elseif proc < 10
            txt = ['.0' int2str(proc)];
        else
            txt = ['.' int2str(proc)];
        end
        text(j, i, txt, 'horizontalalignment', 'center', ...
            'verticalalignment', 'middle', ...
            'fontunits', 'normalized', 'fontsize', params.fontSize/n);
    end
end

end