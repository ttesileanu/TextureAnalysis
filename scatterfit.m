function [c, stats, handles] = scatterfit(x, y, varargin)
% SCATTERFIT Make a scatter plot with a best fit line.
%   [c, stats, handles] = SCATTERFIT(x, y, ...) makes a scatter plot and
%   combines it with a best fit line. This essentially calls smartscatter
%   and drawfitline one after the other, returning c and stats from
%   drawfitline, and combining the handle information from the two
%   functions.
%
%   The options to this function are the options that are allowed for
%   smartscatter and drawfitline.
%
% See also: SMARTSCATTER, DRAWFITLINE.


% parse the optional arguments that scattefit cares about
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('nodraw', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% need to split the options between options meant for smartscatter, and
% options meant for drawfitline; some of them are common
[scatter_opts, fit_opts, other_opts] = splitoptions(varargin, ...
    {'alpha', 'maxPoints'}, ...
    {'conflevel', 'corrtext', 'corrtype', 'intercept', 'legend', ...
        'legendbox', 'legendloc', 'line', 'nodraw', 'showci', 'thinci', 'showfit', ...
        'style', 'refseq', 'r2bootstrap', 'permutation'});
if ~isempty(other_opts)
    % XXX the downside of this is that I need to worry about keeping the
    % list above in-sync with smartscatter and drawfitline
    % an alternative would be to send all unrecognized options to one of
    % them, but which one?
    error([mfilename ':badoption'], ['Unknown option ''' other_opts{1} '''.']);
end

if ~params.nodraw
    handles_scatter = smartscatter(x, y, scatter_opts{:});
    scatter_axes = gca;
    washold = ishold;
    hold on;
else
    handles_scatter = struct;
end

[c, stats, handles_fit] = drawfitline(x, y, fit_opts{:});

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