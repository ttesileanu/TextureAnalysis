function [c, stats, handles] = scatterfit(x, y, varargin)
% SCATTERFIT Make a scatter plot with a best fit line.
%   [c, stats, handles] = SCATTERFIT(x, y) makes a scatter plot and
%   combines it with a best fit line. This essentially calls smartscatter
%   and drawfitline one after the other, returning c and stats from
%   drawfitline, and combining the handle information from the two
%   functions.
%
%   Options
%    'scatteropts'
%       Cell array of options to be passed to smartscatter.
%    'fitopts'
%       Cell array of options to be passed to drawfitline.
%    'nodraw'
%       Don't display any plots, just calculate the statistics from
%       drawfitline.
%
% See also: SMARTSCATTER, DRAWFITLINE.


% parse the optional arguments that scattefit cares about
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('scatteropts', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('fitopts', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('nodraw', false, @(b) isscalar(b) && islogical(b));

% display defaults if asked to
if nargin == 1 && strcmp(x, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

if params.nodraw
    % make sure to pass this to drawfitline, as well
    params.fitopts = [params.fitopts {'nodraw', true}];
    handles_scatter = struct;
else
    handles_scatter = smartscatter(x, y, params.scatteropts{:});
    % keep track of current state so we can return to it after drawfitline
    scatter_axes = gca;
    washold = ishold;
    hold on;
end

% draw the fit line (or get c and stats if 'nodraw' is true)
[c, stats, handles_fit] = drawfitline(x, y, params.fitopts{:});

if ~params.nodraw
    % return to the old state
    axis(scatter_axes);
    if ~washold
        hold off;
    end
end

% collect the handles from both smartscatter and drawfitline
handles = handles_scatter;
fields = fieldnames(handles_fit);
for i = 1:length(fields)
    handles.(fields{i}) = handles_fit.(fields{i});
end

end