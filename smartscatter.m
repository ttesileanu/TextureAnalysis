function handles = smartscatter(x, y, varargin)
% SMARTSCATTER Make a scatter plot with smart subsampling for large
% datasets.
%   SMARTSCATTER(x, y) makes a scatter plot similar to Matlab's scatter
%   command. Unlike that function, SMARTSCATTER automatically flattens the
%   x and y data before plotting, and plots a subset of the data in cases
%   in which the number of points is too large. There are also easy-to-use
%   commands to control transparency (see below).
%
%   SMARTSCATTER uses the same argument format as Matlab's scatter command;
%   use `help scatter` to find out about that. In particular,
%
%   SMARTSCATTER(x, y, s, c) can be used to draw markers of different sizes
%   and colors, just like with Matlab's scatter.
%
%   h = SMARTSCATTER(...) return handles for the scatter objects.
%
%   A few named options specific to SMARTSCATTER are available:
%    'alpha': double
%       Transparency level of markers (0 fully transparent, 1 fully
%       opaque).
%       (default: 0.6)
%     'maxPoints' <n>
%       Maximum number of points to display. If the number of elements in
%       the input vectors is larger than this, a uniform sampling will be
%       done so that only n points are plotted.
%       (default: 5000)
%
%   See also: SCATTER.

% split out named arguments
strMask = cellfun(@(x) ischar(x) && isvector(x), varargin);
strMask(length(strMask):-2:1) = false;
firstNamedArgIdx = find(strMask, 1, 'last');
if ~isempty(firstNamedArgIdx)
    argsPos = varargin(1:firstNamedArgIdx-1);
    argsNamed = varargin(firstNamedArgIdx:end);
else
    argsPos = varargin;
    argsNamed = {};
end
if ~isempty(argsNamed) && length(argsNamed{1}) <= 2
    argsPos = [argsPos argsNamed(1:2)];
    argsNamed = argsNamed(3:end);
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('alpha', 0.6, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('maxPoints', 5000, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(argsNamed{:});
params = parser.Results;

% flatten inputs
x = x(:);
y = y(:);

n = length(x);
if length(y) ~= n
    error([mfilename ':badsz'], 'Numbers of elements of x and y must match.');
end

% subsample if necessary
if n > params.maxPoints
    idxs = round(linspace(1, n, params.maxPoints));
    x = x(idxs);
    y = y(idxs);
    
    for i = 1:length(argsPos)
        % XXX this is very error-prone!
        if isvector(argsPos{i}) && length(argsPos{i}) == n
            argsPos{i} = argsPos{i}(idxs);
        elseif ismatrix(argsPos{i}) && size(argsPos{i}, 1) == n
            argsPos{i} = argsPos{i}(idxs, :);
        end
    end
end

% figure out scatter options
if params.alpha < 1
    opts = {'markerfacealpha', params.alpha, 'markeredgealpha', params.alpha};
else
    opts = {};
end
opts = [opts structToCell(parser.Unmatched)];

% plot!
handles = scatter(x, y, argsPos{:}, opts{:});

end