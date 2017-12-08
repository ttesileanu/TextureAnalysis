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

% parse
parser.parse(varargin{:});
params = parser.Results;

if params.nodraw
    % make sure to pass this to drawfitline, as well
    params.fitopts = [params.fitopts {'nodraw', true}];
end

if ~params.nodraw
    handles_scatter = smartscatter(x, y, params.scatteropts{:});
    scatter_axes = gca;
    washold = ishold;
    hold on;
else
    handles_scatter = struct;
end

[c, stats, handles_fit] = drawfitline(x, y, params.fitopts{:});

if ~params.nodraw
    % the axes should be given by the scatter plot, not by the fit line
    axis(scatter_axes);
    if ~washold
        hold off;
    end
end

handles = handles_scatter;
fields = fieldnames(handles_fit);
for i = 1:length(fields)
    handles.(fields{i}) = handles_fit.(fields{i});
end

end