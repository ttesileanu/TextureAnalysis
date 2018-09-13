function h = textTexGroup(x, y, group, varargin)
% textTexGroup Write a texture group symbol in an axis.
%   textTexGroup(x, y, group) displays the name of the `group` at location
%   `(x, y)` in the current axes in the format used in the paper.
%
%   Note that although the position of the text is by default given in axes
%   coordinates, it will always stay fixed *in pixel coordinates*, so any
%   resizing of the figure will change the axis position of the text.
%
%   h = textTexGroup(...) returns handles to the text objects that were
%   created.
%
%   Options (case insensitive to match Matlab's TEXT):
%    'FontName'
%       Name of font to use.
%    'FontSize'
%       Font size for the top-level symbol.
%    'Units':
%       Units in which the `x` and `y` arguments are given (default is data
%       units).
%    'HorizontalAlignment'
%       Set how the text is aligned horizontally. This can be 'left',
%       'center', or 'right'.
%    'VerticalAlignment'
%       Set how the text is aligned vertically. This can be 'middle',
%       'top', or 'bottom'. Note that 'middle' refers to the *symbol only*,
%       while 'top' and 'bottom' refer to the whole group name.
%    'symbol'
%       Symbol to use for the texture group. If not given, it uses '\gamma'
%       for first-order group, '\beta' for second order, '\theta' for third
%       order, and '\alpha' for fourth order.
%    'displayZeros'
%       If `true`, displays even coefficients that are 0.
%    'postfix'
%       Text to write after the group name and coefficients.
%    'coeffToStr'
%       Custom function to use to convert coefficients to strings.
%    'subscriptScaling'
%       Size of subscript font as a fraction of the symbol font.
%    'subscriptSpacing'
%       Extra spacing added between subscripts, in units of subscript
%       extents. Set to a negative value to diminish spacing.
%    'postfixSpacing'
%       Spacing added before postfix, in units of subscript extents.
%    'groupSpacing'
%       Spacing added between groups, in units of subscript extents.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = false;
parser.FunctionName = mfilename;

parser.addParameter('FontName', '', @(s) ischar(s) && isvector(s));
parser.addParameter('FontSize', [], @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('Units', '', @(s) ischar(s) && isvector(s));
parser.addParameter('HorizontalAlignment', 'left', @(s) ismember(s, {'left', 'center', 'right'}));
parser.addParameter('VerticalAlignment', 'middle', @(s) ismember(s, {'middle', 'top', 'bottom'}));
parser.addParameter('symbol', '', @(s) ischar(s) && isvector(s));
parser.addParameter('displayZeros', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('postfix', '', @(s) ischar(s) && isvector(s));
parser.addParameter('coeffToStr', @int2str);
parser.addParameter('subscriptScaling', 0.8, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('subscriptSpacing', -0.4, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('postfixSpacing', 0.1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('groupSpacing', 0.2, @(x) isnumeric(x) && isscalar(x));

% show defaults if requested
if nargin == 1 && strcmp(x, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% check whether this is a mixed group
nGroups = 1 + sum(group == ';');

h = [];
if nGroups > 1
    groups = strsplit(group, ';');
else
    groups = {group};
end

units = params.Units;

for idx = 1:length(groups)
    group = groups{idx};
    % parse the group name
    
    % get rid of useless whitespace
    group = strtrim(group);
    
    % split out direction specification, if it exsists
    bracketIdx = find(group == '[', 1);
    if ~isempty(bracketIdx)
        bracket = group(bracketIdx:end);
        group = group(1:bracketIdx-1);
    else
        bracket = '';
    end
    
    % split group type and coefficients
    parts = strsplit(group, '_');
    type = parts{1};
    coefficients = cellfun(@str2double, parts(2:end));
    
    % find group order, identify letter to write
    if isempty(params.symbol)
        order = length(type);
        symbolList = {'\gamma', '\beta', '\theta', '\alpha'};
        symbol = symbolList{order};
    else
        symbol = params.symbol;
    end
    
    % print the symbol name
    textOptions = {};
    if ~isempty(params.FontName)
        textOptions = [textOptions {'FontName', params.FontName}]; %#ok<AGROW>
    end
    if ~isempty(params.FontSize)
        textOptions = [textOptions {'FontSize', params.FontSize}]; %#ok<AGROW>
    end
    if ~isempty(units)
        textOptions = [textOptions {'Units', units}]; %#ok<AGROW>
    end
    hSymbol = text(x, y, symbol, textOptions{:});
    
    % figure out size and spacing for subscripts
    subscriptSize = hSymbol.FontSize*params.subscriptScaling;
    
    tmp = uicontrol('style', 'text', 'string', 'X', ...
        'fontname', hSymbol.FontName, 'fontsize', subscriptSize, ...
        'visible', false);
    % note that the control is on the figure, so we can't get extents in axis
    % coordinates
    tmp.Units = 'pixels';
    subscriptExtents = tmp.Extent(3:4);
    % the extents obtained like this are way too big...
    subscriptExtents = subscriptExtents*(1 + params.subscriptSpacing);
    
    % specify position of subscripts in pixels
    hSymbol.Units = 'pixels';
    symbolEdge = hSymbol.Position(1:2) + [hSymbol.Extent(3) 0];
    
    % figure out coefficients for all positions (A, B, C, and D)
    allCoeffs = zeros(1, 4);
    allCoeffs(type - 'A' + 1) = coefficients;
    
    % draw the coefficients
    largestX = 0;
    hCoeffs = [];
    for k = 1:4
        i = floor((k - 1)/2);
        j = mod(k - 1, 2);
        if allCoeffs(k) > 0 || params.displayZeros
            crtH = text(symbolEdge(1) + j*subscriptExtents(1), symbolEdge(2) - i*subscriptExtents(2), ...
                params.coeffToStr(allCoeffs(k)), ...
                'fontname', hSymbol.FontName, 'fontsize', subscriptSize, ...
                'verticalalignment', 'top', 'Units', 'pixels');
            largestX = max(largestX, j*subscriptExtents(1));
            hCoeffs = [hCoeffs crtH]; %#ok<AGROW>
        end
    end
    
    % draw the bracket and/or postfix, if any
    if idx == length(groups)
        postfix = strcat(bracket, params.postfix);
    else
        % if this isn't the last group, draw a semicolon as a postfix
        postfix = strcat(bracket, ';');
    end
    if ~isempty(postfix)
        bracketPosition = symbolEdge + ...
            [largestX + subscriptExtents(1)*(1 + params.postfixSpacing), 0];
        crtH = text(bracketPosition(1), bracketPosition(2), postfix, 'fontname', ...
            hSymbol.FontName, 'fontsize', hSymbol.FontSize, 'Units', 'pixels');
        hCoeffs = [hCoeffs crtH]; %#ok<AGROW>
    end
    
    h = [h hSymbol hCoeffs]; %#ok<AGROW>
    
    % move to next group
    x = h(end).Extent(1) + h(end).Extent(3) + subscriptExtents(1)*params.groupSpacing;
    y = h(end).Position(2);
    
    % subsequent units are always pixels
    units = 'pixels';
end

if ~strcmp(params.HorizontalAlignment, 'left') || ~strcmp(params.VerticalAlignment, 'middle')
    % find the full extents of the group name
    allExtent = cell2mat({h.Extent}');
    
    % convert to [min_x, min_y, max_x, max_y]
    actualExtents = [allExtent(:, 1:2) allExtent(:, 1:2) + allExtent(:, 3:4)];
    
    fullExtents = [min(actualExtents(:, 1:2), [], 1) max(actualExtents(:, 3:4), [], 1)];
    
    dimensions = fullExtents(3:4) - fullExtents(1:2);
    
    shift = [0 0];
    switch params.HorizontalAlignment
        case 'center'
            shift(1) = -dimensions(1) / 2;
        case 'right'
            shift(1) = -dimensions(1);
    end
    
    % get original y value in pixels
    h(1).Units = 'pixels';
    originalY = h(1).Position(2);
    
    switch params.VerticalAlignment
        case 'top'
            shift(2) = originalY - fullExtents(4);
        case 'bottom'
            shift(2) = originalY - fullExtents(2);
    end

    for i = 1:length(h)
        h(i).Units = 'pixels'; %#ok<AGROW>
        h(i).Position(1:2) = h(i).Position(1:2) + shift; %#ok<AGROW>
    end
end

end
