function preparegraph(varargin)
% PREPAREGRAPH Set the figure's paper size so that it can be properly
% exported to PDF.
%
%   PREPAREGRAPH(handle) sets the figure's paper size so that it can be
%   properly exported to PDF.
%
%   PREPAREGRAPH acts on the current figure.
%
%   Options:
%    'edge' <f>
%       Fraction of width/height to use as edge.
%       (default: 0.01)

% parse the optional arguments that scattefit cares about
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('handle', gcf);
parser.addParameter('edge', 0.01, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

oldUnits = get(params.handle, 'Units');
units = get(params.handle, 'PaperUnits');
set(params.handle, 'Units', units);
screenpos = get(params.handle, 'Position');
w = screenpos(3);
h = screenpos(4);
set(params.handle, 'Units', oldUnits);

set(params.handle, ...
    'PaperPosition', [params.edge*w, params.edge*h, w, h], ...
    'PaperSize', [w*(1 + 2*params.edge), h*(1 + 2*params.edge)] ...
);

end