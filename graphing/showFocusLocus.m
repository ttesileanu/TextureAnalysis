function showFocusLocus(res, varargin)
% showFocusLocus Show the locus of patches identified as in-focus.
%   showFocusLocus(res) makes scatter plots for all the 10 dimensions of
%   statistics, showing the locus of patches that were identified as
%   in-focus.
%
%   The function ignores spurious patches for which all of the statistics
%   have been set to 1.
%
%   Options:
%    'fSel' <v>
%       Binary mask showing which patches to be considered in-focus.
%       (default: inferred from res.focus.clusterIds and res.focus.focusCluster)
%    'alpha' <x>
%       Transparency value for scatter plot.
%    'size' <x>
%       Size to use for scatter plot.

% parse the optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('fSel', [], @(v) isvector(v) && islogical(v));
parser.addParameter('alpha', [], @(x) isscalar(x) && isnumeric(x) && isreal(x));
parser.addParameter('size', 10, @(x) isscalar(x) && isnumeric(x) && isreal(x));

parser.parse(varargin{:});
params = parser.Results;

%figure;
set(gcf, 'position', [50 50 961 600], 'paperposition', [0.25 0.25 9.61 6]);

% get rid of spurious patches
mask = (res.ev(:, 1) ~= 1);

focus_comp = res.focus.focusCluster;
if isempty(params.fSel)
    mask_focus_comp = (res.focus.clusterIds == focus_comp);
else
    mask_focus_comp = params.fSel;
end
mask_focus = (mask(:) & mask_focus_comp(:));
mask_blur = (mask(:) & ~mask_focus_comp(:));

c_labels = {'\gamma', '\beta_|', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{|-}', '\theta_{-|}', '\theta_{\_|}', '\theta_{|\_}', '\alpha'};

order = [1, 10, 2, 3, 4, 5, 6, 7, 8, 9];

if isempty(params.alpha)
    opts = {};
else
    opts = {'alpha', params.alpha};
end

for i0 = 1:2:10
    i1 = order(i0);
    i2 = order(i0+1);
    
    subplot(2, 3, (i0+1)/2);
    
    smartscatter(res.ev(mask_blur, i1), res.ev(mask_blur, i2), params.size, ...
        'filled', opts{:});
    xl1 = xlim;
    yl1 = ylim;
    hold on;
    
    smartscatter(res.ev(mask_focus, i1), res.ev(mask_focus, i2), params.size, ...
        'filled', opts{:});
    xl2 = xlim;
    yl2 = ylim;
    
    xl = [min(xl1(1), xl2(1)) max(xl1(2), xl2(2))];
    yl = [min(yl1(1), yl2(1)) max(yl1(2), yl2(2))];
    
    xlim(xl);
    ylim(yl);
    
    legend('blurred', 'in-focus');
    
    xlabel(['comp. ' int2str(i1) ': ' c_labels{i1}]);
    ylabel(['comp. ' int2str(i2) ': ' c_labels{i2}]);
    
%    legend({'in-focus', 'blurred'});
end

end