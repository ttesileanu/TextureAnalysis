function beautifygraph(varargin)
% BEAUTIFYGRAPH Improve the quality of the graph.
%   BEAUTIFYGRAPH improves the quality of the current active figure by
%   changing font styles and sizes and tick directions.
%
%   BEAUTIFYGRAPH(haxes) is applied to the axes identified by the given
%   handle.
%
%   Options:
%    'changefonts' <b>
%       True to change fonts.
%       (default: false)
%    'fontscale' <x>
%       A scale factor for all the font sizes.
%       (default: 1)

% This function is based on the suggestions given at
% blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs

% Tiberiu Tesileanu (2013-2014)

% parse the optional arguments that scattefit cares about
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('axes', gca);
parser.addParamValue('fontscale', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('changefonts', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

htitle = get(params.axes, 'title');
hxlabel = get(params.axes, 'xlabel');
hylabel = get(params.axes, 'ylabel');

if params.changefonts
    set(params.axes, 'fontname', 'Helvetica');
    set([htitle, hxlabel, hylabel], 'fontname', 'AvantGarde');
end
set([hxlabel, hylabel], 'fontsize', 12*params.fontscale);
set(htitle, 'fontsize', 14*params.fontscale, 'fontweight', 'bold');
% XXX can you rescale font for tick labels?

set(params.axes, ...
    'box',          'off', ...
    'tickdir',      'out', ...
    'xminortick',   'on', ...
    'yminortick',   'on', ...
    'linewidth',    1);

end