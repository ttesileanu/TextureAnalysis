function chisqp = multihist(varargin)
% MULTIHIST Overlap several histograms.
%   MULTIHIST(ys) plots one or more histograms. The input can be a single
%   vector, if one histogram is to be plotted, or a cell array of vectors,
%   one per histogram. The bin centers will be determined automatically. By
%   default this uses 10 bins.
%
%   MULTIHIST(ys, n) specifies the number of bins to use.
%
%   MULTIHIST(ys, x) plots one or more histograms, using the provided bin
%   centers.
%
%   chisqp = MULTIHIST(...) returns a matrix of chi-squared test p-values
%   for all the pairs of histograms (after scaling). The (i, j) entry is
%   the value of the test assuming j are the expected values and i are the
%   observed ones.
%
%   Options:
%    'colors' <v>
%       Set the colors used by the histograms; one option is to use a
%       string of characters following Matlab's plot color options.
%       Alternatively, this can be a k x 3 real vector, where each
%       row is an RGB-specified color. In both cases, the colors are used
%       cyclically if there are more histograms than colors.
%       (default: 'brg')
%    'mode' <s>
%       This can be
%         'overlap':    draw the histograms on top of each other
%         'stack':      make a stacked histogram
%       (default: 'overlap')
%    'offset' <x/v>
%       How much to offset the bars, as a fraction of a bar width. This can
%       be a scalar to offset all histograms by the same amount or, more
%       usefully, a vector, to give different offsets to different
%       histograms. This does not work when 'mode' is 'stack'.
%       (default: 0)
%    'scaling' <s/v>
%       Choose whether to scale the histograms before displaying; options
%       are a vector of scale factors for each of the histograms, or one of
%       the following:
%           'none' -- no scaling is performed
%           'max'  -- the maximum heights of the two histograms are aligned
%           'mean' -- the mean heights of the two histograms are aligned.
%       (default: 'none')
%    'width' <w/v>
%       Sets the width of the histogram bars, as used by Matlab's bar
%       command. If the 'mode' isn't 'stack', this can be a vector, to
%       impose independent settings for each histogram.
%       (default: 1)

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addRequired('ys', @(v) isvector(v) && (iscell(v) || isnumeric(v)));
parser.addOptional('x', [], @(x) isvector(x) && isnumeric(x));

parser.addParamValue('colors', 'brg', @(s) (ischar(s) && isvector(s)) || (isnumeric(s) && ismatrix(s) && size(s, 2) == 3));
parser.addParamValue('mode', 'overlap', @(s) ismember(s, {'overlap', 'stack'}));
parser.addParamValue('offset', 0, @(v) isnumeric(v) && isvector(v));
parser.addParamValue('scaling', 'none', @(v) isvector(v) && (isnumeric(v) || (ischar(v) && ismember(v, {'none', 'max', 'mean'}))));
parser.addParamValue('width', 1, @(v) isnumeric(v) && isvector(v));

% parse
parser.parse(varargin{:});
params = parser.Results;

x = params.x;
ys = params.ys;
params.nbins = 10;

% check that the data is numeric
if ~isempty(x)
    if isscalar(x)
        if x ~= fix(x)
            error([mfilename ':nonint'], 'The number of bins should be an integer.');
        end
        params.nbins = x;
        x = [];
    else
        if (~isnumeric(x) && ~islogical(x)) || ~isvector(x)
            error([mfilename ':badx'], 'Bin center data should be a numeric vector.');
        end
    end
end
if ~all(cellfun(@(a) (isnumeric(a) || islogical(a)) && isvector(a), ys))
    error([mfilename ':badys'], 'All the data should be numeric vectors.');
end

% handle scalar vs. vector width and offset
n = length(ys);
if isscalar(params.width)
    params.width = repmat(params.width, n, 1);
end
if isscalar(params.offset)
    params.offset = repmat(params.offset, n, 1);
end

if isempty(x)
    % decide bin centers
    minval = min(cellfun(@min, ys));
    maxval = max(cellfun(@max, ys));
    step = (maxval - minval) / (params.nbins - 1);
    
    x = minval:step:maxval;
end

% create the histograms
histos = cell(1, n);
for i = 1:n
    histos{i} = flatten(hist(ys{i}, x));
end

% scale if we were asked to
if ischar(params.scaling)
    switch params.scaling
        case 'max'
            maxima = cellfun(@max, histos);
            overallmax = max(maxima);
            for i = 1:n
                histos{i} = histos{i} * overallmax / maxima(i);
            end
        case 'mean'
            means = cellfun(@mean, histos);
            overallmean = mean(means);
            for i = 1:n
                histos{i} = histos{i} * overallmean / means(i);
            end
        case 'none'
    end
else
    if ~isvector(params.scaling) || ~isnumeric(params.scaling) || ~isreal(params.scaling)
        error([mfilename ':badscalingvec'], 'Scaling must be either one of the string options or a vector of reals.');
    end
    for i = 1:n
        histos{i} = histos{i}*params.scaling(i);
    end
end

if nargout >= 1
    chisqp = zeros(n, n);
    for i = 1:n
        for j = 1:n
            if j == i
                chisqp(i, j) = 1;
            else
                chisqp(i, j) = chi2gof2([histos{i} histos{j}]);
            end
        end
    end
end

% draw the histograms
washold = ishold;
hold on;
if ischar(params.colors)
    ncols = length(params.colors);
else
    ncols = size(params.colors, 1);
end
set(gca, 'fontsize', 12);

if strcmp(params.mode, 'stack')
    params.colors = torgb(params.colors);
    mycolormap = [...
        repmat(params.colors, floor(n/ncols), 1) ; ...
        params.colors(1:mod(n, ncols), :) ...
    ];
    
    % can only use one width here!
    bar(x, cell2mat(histos), params.width(1), 'stacked');
    colormap(mycolormap);
else
    xw = x(2) - x(1);
    for i = 1:n
        if ischar(params.colors)
            color = params.colors(mod(i - 1, ncols) + 1);
        else
            color = params.colors(mod(i - 1, ncols) + 1, :);
        end
        % XXX aargh! stupid Matlab 2014b outputs a weird object when using
        % the non-'hist' version of bar, and this weird object does not
        % support transparency...
        tmp = bar(x + params.offset(i)*xw, histos{i}, params.width(i), 'hist');
        set(tmp, 'facecolor', color);
    end
end

if ~washold
    hold off;
end

end