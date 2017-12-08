function [c, stats, handles] = drawfitline(x, y, varargin)
% DRAWFITLINE Draw a best fit line.
%   DRAWFITLINE(x, y) calculates a best fit line for the given points and
%   draws it in the current axes. The function automatically flattens the
%   data and ignores NaNs.
%
%   c = DRAWFITLINE(...) returns the correlation coefficient between x and y.
%
%   [c, stats] = DRAWFITLINE(...) also returns a structure with other
%   statistical information, including the best fit parameters a and b
%   (such that y = a*x + b), the p-value for the null hypothesis that the
%   slope is zero, and a matrix of confidence intervals for a and b.
%
%   [c, stats, handles] = DRAWFITLINE(...) returns a structure of handles
%   for the fit line and/or the confidence intervals.
%
%   Options:
%     'conflevel' <x>
%       Confidence level to use for error estimates.
%       (default: 0.95)
%     'corrtext' <s>
%       Text to write when displaying the correlation coefficient.
%       (default: 'c = ' or [corrtype ' c = '])
%     'corrtype' <s>
%       Type of correlation coefficient to return. This should be one of the
%       options supported by Matlab's corr function.
%       (default: Matlab's default)
%     'intercept' <x>
%       Fix the intercept of the fit line to x.
%       (default: treat the intercept as a fitting parameter)
%     'legend' <s/b>
%       What to display in the legend, or false to not display a legend.
%       This should be a string in which each character shows a piece of
%       information to display. These are the allowed options:
%         'c':  correlation coefficient
%         'f':  formula for the fit line
%         'i':  confidence intervals for fit parameters
%         'p':  p-value for null hypothesis that slope is zero
%         's':  standard deviation of residuals
%       (default: 'fc')
%     'legendbox' <b>
%       Whether the legend should have a bounding box.
%       (default: false)
%     'legendloc' <x>
%       Location for the legend. This can be 'north', 'northeast', 'east',
%       'southeast', 'south', 'southwest', 'west', or 'northwest'.
%       (default: 'north')
%     'line' [a, b]
%       Force the line to be drawn with slope a and intercept b. This does
%       not interfere with the return values of this function.
%       (default: draw the best fit line)
%     'nodraw' <b>
%       If true, do not display the line; instead, just return the fit
%       coefficients and error estimates.
%       (default: false)
%     'permutation' <n>
%       If nonzero, sets the number of samples used in a permutation test
%       to establish the (two-sided) p-value for the null hypothesis that
%       there is no correlation between the two data series.
%     'r2bootstrap' <n>
%       If nonzero, sets the number of bootstrap samples used to estimate a
%       confidence interval for the correlation coefficient.
%     'showci' <b>
%       Whether to draw the confidence interval for prediction on the plot.
%       (default: true)
%     'thinci' <b>
%       When 'showci' is true, whether to show a thin confidence interval,
%       that does not contain the interval for the residuals.
%       (default: true)
%     'skipinf' <b>
%       If true, skip values that are infinite or NaN in either vector.
%       (default: true)
%     'showfit' <b>
%       Whether to show the fit line or not. Note that the confidence
%       interval can be displayed independently of the line.
%       (default: true)
%     'style' <c>
%       Arguments to pass to plot for fit line.
%       (default: {'k', 'linewidth', 2})
%
% See also: CORR.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('conflevel', 0.95, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('corrtext', '', @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParameter('corrtype', '', @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParameter('intercept', [], @(x) isscalar(x) && isnumeric(x));
parser.addParameter('legend', 'fc', @(b) (islogical(b) && isscalar(b) && ~b) || (ischar(b) && isvector(b)));
parser.addParameter('legendbox', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('legendloc', 'north', @(s) ismember(s, {'north', 'northeast', 'east', 'southeast', 'south', 'southwest', 'west', 'northwest'}));
parser.addParameter('line', [], @(v) isvector(v) && isnumeric(v) && length(v) == 2);
parser.addParameter('nodraw', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('showci', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('thinci', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('skipinf', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('showfit', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('style', {'k', 'linewidth', 2}, @(c) iscell(c));

parser.addParameter('r2bootstrap', 0, @(n) isscalar(n) && isnumeric(n));
parser.addParameter('permutation', 0, @(n) isscalar(n) && isnumeric(n));

% parse
parser.parse(varargin{:});
params = parser.Results;

x = x(:);
y = y(:);

if length(x) ~= length(y)
    error([mfilename ':badsizes'], 'The sizes of x and y do not match.');
end

if params.skipinf
    mask = isfinite(x) & isfinite(y);
end
x = x(mask);
y = y(mask);

% calculate the correlation coefficient and the p-value
if isempty(params.corrtype)
    corrargs = {};
else
    corrargs = {'type', params.corrtype};
end
[c, p] = corr(x, y, corrargs{:});

% calculate the best fit line
if isempty(params.intercept)
    [coeffs, fitstruc] = polyfit(x, y, 1);
    a = coeffs(1);
    b = coeffs(2);
    stats = fitstruc;
else
    a = x\(y - params.intercept);
    b = params.intercept;
end

stats.a = a;
stats.b = b;
stats.p = p;
stats.c = c;

% handle bootstrap estimate for confidence interval of correlation coeff.
if params.r2bootstrap > 0
    bs_x = cell(params.r2bootstrap, 1);
    bs_y = cell(params.r2bootstrap, 1);
    bs_c = zeros(params.r2bootstrap, 1);
    bs_p = zeros(params.r2bootstrap, 1);
    
    n = length(x);
    
    for i = 1:params.r2bootstrap
        choices = randi([1 n], n, 1);
        bs_x{i} = x(choices);
        bs_y{i} = y(choices);
        [bs_c(i), bs_p(i)] = corr(bs_x{i}, bs_y{i}, corrargs{:});
    end
    
    stats.bootstrap_c_all = bs_c;
    stats.bootstrap_c_mean = mean(bs_c);
    stats.bootstrap_c_std = std(bs_c);
    
    stats.bootstrap_p_all = bs_p;
    stats.bootstrap_p_mean = mean(bs_p);
    stats.bootstrap_p_std = std(bs_p);
end

% handle permutation test estimate for p value
if params.permutation > 0
    perm_y = cell(params.permutation, 1);
    perm_c = zeros(params.permutation, 1);
    perm_p = zeros(params.permutation, 1);
    
    n = length(x);
    
    for i = 1:params.permutation
        perm_y{i} = y(randperm(n, n));
        [perm_c(i), perm_p(i)] = corr(x, perm_y{i}, corrargs{:});
    end
    
    stats.permutation_p = mean(abs(perm_c) > abs(c));
    stats.permutation_c_all = perm_c;
end

% calculate residuals and confidence intervals
% XXX this might not be quite right when the intercept is fixed
residuals = y - a*x - b;
varresiduals = sum(residuals.^2) / max(1, length(x)-2);
sumx2centered = sum((x - mean(x)).^2);
sa = sqrt(varresiduals / sumx2centered);
sb = sa*sqrt(mean(x.^2));
qstar = tinv((1 + params.conflevel) / 2, length(x) - 2);
erra = qstar*sa;
errb = qstar*sb;

cis = zeros(2, 2);

cis(1, 1) = a - erra;
cis(2, 1) = a + erra;
if isempty(params.intercept)
    cis(1, 2) = b - errb;
    cis(2, 2) = b + errb;
else
    cis(:, 2) = b;
end

handles = struct;
if ~params.nodraw
    dummyx = [min(x) max(x)];
    if isempty(params.line)
        desc = 'fit';
        pa = a;
        pb = b;
    else
        desc = 'guide';
        pa = params.line(1);
        pb = params.line(2);
    end
    
    dummyy = pa*dummyx + pb;
    if params.showfit
        handles.fit = plot(dummyx, dummyy, params.style{:});
    end
        
    if params.showci && strcmp(desc, 'fit')
        % draw the confidence interval
        avgx = mean(x);
        dummyx = linspace(min(x), max(x), 1000);
        dummyy = pa*dummyx + pb;
        % XXX I'm a little fuzzy on where exactly the 1/length(x) term
        % comes from, but it seems to show up in people's formulas
        intervals = qstar*sqrt(varresiduals*(~params.thinci + 1/length(x) + (dummyx - avgx).^2 / sumx2centered));
        ci_x = [dummyx fliplr(dummyx)];
        ci_y = [(dummyy-intervals) fliplr(dummyy+intervals)];
        handles.ci = fill(ci_x, ci_y, [1 0.8 0.8], 'edgecolor', 'none');
        uistack(handles.ci, 'bottom');
    end
    
    if ischar(params.legend)
        todisp = {};
        for i = 1:length(params.legend)
            switch params.legend(i)
                case 'c'
                    if params.r2bootstrap == 0
                        cvalue = num2str(c, 3);
                    else
                        cvalue = [num2str(stats.bootstrap_c_mean, 2) '+-' ...
                            num2str(stats.bootstrap_c_std, 2)];
                    end
                    if isempty(params.corrtext)
                        ctype = lower(params.corrtype);
                        if ~isempty(ctype)
                            ctype(1) = upper(ctype(1));
                            ctype = [ctype ' ']; %#ok<AGROW>
                        end
                        crttext = [ctype 'c = ' cvalue];
                    else
                        crttext = [params.corrtext cvalue];
                    end
                case 'f'
                    crttext = [desc ' y = ' num2str(pa, 3) ' x '];
                    if pb >= 0
                        crttext = [crttext '+']; %#ok<AGROW>
                    else
                        crttext = [crttext '-']; %#ok<AGROW>
                    end
                    crttext = [crttext ' ' num2str(abs(pb), 3)]; %#ok<AGROW>
                case 'i'
                    crttext = {'95% CI:', ...
                        ['    a=' num2str(cis(1, 1), 3) ' -- ' num2str(cis(2, 1), 3)], ...
                        ['    b=' num2str(cis(1, 2), 3) ' -- ' num2str(cis(2, 2), 3)] ...
                    };
                case 'p'
                    crttext = ['p = ' num2str(p, 2)];
                case 's'
                    crttext = ['std = ' num2str(sqrt(varresiduals), 2)];
                otherwise
                    error([mfilename ':badlt'], ['Unrecognized legend character ''' params.legend(i) '''.']);
            end
            todisp = [todisp crttext]; %#ok<AGROW>
        end
        
        len = max(cellfun(@length, todisp)) + 1;
        descstr = repmat(blanks(len), length(todisp), 1);
        for i = 1:length(todisp)
            descstr(i, 1:1+length(todisp{i})) = [' ' todisp{i}];
        end
        
        tmp = axis;
        switch params.legendloc
            case 'north'
                txtpos = [mean(tmp(1:2)) tmp(4)];
                txtprops = {'horizontalalignment', 'center', ...
                    'verticalalignment', 'top'};
            case 'northeast'
                txtpos = [tmp(2) tmp(4)];
                txtprops = {'horizontalalignment', 'right', ...
                    'verticalalignment', 'top'};
            case 'northwest'
                txtpos = [tmp(1) tmp(4)];
                txtprops = {'horizontalalignment', 'left', ...
                    'verticalalignment', 'top'};
            case 'south'
                txtpos = [mean(tmp(1:2)) tmp(3)];
                txtprops = {'horizontalalignment', 'center', ...
                    'verticalalignment', 'bottom'};
            case 'southeast'
                txtpos = [tmp(2) tmp(3)];
                txtprops = {'horizontalalignment', 'right', ...
                    'verticalalignment', 'bottom'};
            case 'southwest'
                txtpos = [tmp(1) tmp(3)];
                txtprops = {'horizontalalignment', 'left', ...
                    'verticalalignment', 'bottom'};
            case 'east'
                txtpos = [tmp(2) mean(tmp(3:4))];
                txtprops = {'horizontalalignment', 'right', ...
                    'verticalalignment', 'middle'};
            case 'west'
                txtpos = [tmp(1) mean(tmp(3:4))];
                txtprops = {'horizontalalignment', 'left', ...
                    'verticalalignment', 'middle'};
            otherwise
                error([mfilename ':badloc'], ['Unsupported legend location: ' params.legend_location '.']);
        end
        
        if params.legendbox
            txtprops = [txtprops {'edgecolor', [0 0 0]}];
        end
        
        if ~isempty(descstr)
            handles.legend = text(txtpos(1), txtpos(2), descstr, txtprops{:});
        end        
    end    
end

stats.ci = cis;
stats.residuals = residuals;

end