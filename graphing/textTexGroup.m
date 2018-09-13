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
%    'FontName':
%       Name of font to use.
%    'FontSize':
%       Font size for the top-level symbol.
%    'Units':
%       Units in which the `x` and `y` arguments are given (default is data
%       units).
%    'symbol':
%       Symbol to use for the texture group. If not given, it uses '\gamma'
%       for first-order group, '\beta' for second order, '\theta' for third
%       order, and '\alpha' for fourth order.
%    'displayZeros':
%       If `true`, displays even coefficients that are 0.
%    'postfix':
%       Text to write after the group name and coefficients.
%    'coeffToStr':
%       Custom function to use to convert coefficients to strings.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = false;
parser.FunctionName = mfilename;

parser.addParameter('FontName', '', @(s) ischar(s) && isvector(s));
parser.addParameter('FontSize', [], @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('Units', '', @(s) ischar(s) && isvector(s));
parser.addParameter('symbol', '', @(s) ischar(s) && isvector(s));
parser.addParameter('displayZeros', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('postfix', '', @(s) ischar(s) && isvector(s));
parser.addParameter('coeffToStr', @int2str);

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
    subscriptSize = hSymbol.FontSize*0.8;
    
    tmp = uicontrol('style', 'text', 'string', 'X', ...
        'fontname', hSymbol.FontName, 'fontsize', subscriptSize, ...
        'visible', false);
    % note that the control is on the figure, so we can't get extents in axis
    % coordinates
    tmp.Units = 'pixels';
    subscriptExtents = tmp.Extent(3:4);
    % the extents obtained like this are way too big...
    subscriptExtents = subscriptExtents*0.6;
    
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
            h = text(symbolEdge(1) + j*subscriptExtents(1), symbolEdge(2) - i*subscriptExtents(2), ...
                params.coeffToStr(allCoeffs(k)), ...
                'fontname', hSymbol.FontName, 'fontsize', subscriptSize, ...
                'verticalalignment', 'top', 'Units', 'pixels');
            largestX = max(largestX, j*subscriptExtents(1));
            hCoeffs = [hCoeffs h]; %#ok<AGROW>
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
        bracketPosition = symbolEdge + [largestX + subscriptExtents(1)*1.1, 0];
        h = text(bracketPosition(1), bracketPosition(2), postfix, 'fontname', ...
            hSymbol.FontName, 'fontsize', hSymbol.FontSize, 'Units', 'pixels');
        hCoeffs = [hCoeffs h]; %#ok<AGROW>
    end
    
    h = [h hSymbol hCoeffs]; %#ok<AGROW>
    
    % move to next group
    x = h(end).Extent(1) + h(end).Extent(3) + subscriptSize*0.2;
    y = h(end).Position(2);
    
    % subsequent units are always pixels
    units = 'pixels';
end

end
