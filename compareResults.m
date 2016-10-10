function compareResults(res1, res2, varargin)
% compareResults Compare the distributions of texture statistics between
% two image sets.
%   compareResults(res1, res2) makes scatter plots comparing the
%   distributions of texture statistics in the two sets of results.
%
%   The function ignores spurious patches for which all of the statistics
%   have been set to 1.
%
%   Options:
%    'labels' <c>
%       Cell array of labels to use for the legend.
%       (default: no legend)
%    'fields' {<s1>, <s2>}
%       Fields to use from each result structure.
%       (default: {'ev', 'ev'})

% parse the optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('labels', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('fields', {'ev', 'ev'}, @(c) iscell(c) && numel(c) == 2);

parser.parse(varargin{:});
params = parser.Results;

%figure;
set(gcf, 'position', [50 50 961 600], 'paperposition', [0.25 0.25 9.61 6]);

% get rid of spurious patches
ev1 = res1.(params.fields{1});
ev2 = res2.(params.fields{2});
mask1 = (ev1(:, 1) ~= 1);
mask2 = (ev2(:, 1) ~= 1);

c_labels = {'\gamma', '\beta_|', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{|-}', '\theta_{-|}', '\theta_{\_|}', '\theta_{|\_}', '\alpha'};

order = [1, 10, 2, 3, 4, 5, 6, 7, 8, 9];

for i0 = 1:2:10
    i1 = order(i0);
    i2 = order(i0+1);
    
    subplot(2, 3, (i0+1)/2);
    smartscatter(ev1(mask1, i1), ev1(mask1, i2), 10, 'b');
    xl1 = xlim;
    yl1 = ylim;

    hold on;
    smartscatter(ev2(mask2, i1), ev2(mask2, i2), 10, 'r');
    xl2 = xlim;
    yl2 = ylim;
    
    xl = [min(xl1(1), xl2(1)) max(xl1(2), xl2(2))];
    yl = [min(yl1(1), yl2(1)) max(yl1(2), yl2(2))];
    
    xlim(xl);
    ylim(yl);
    
    xlabel(['comp. ' int2str(i1) ': ' c_labels{i1}]);
    ylabel(['comp. ' int2str(i2) ': ' c_labels{i2}]);
    
    if ~isempty(params.labels)
        legend(params.labels{:});
    end
    
%    xlim([-1 1]);
%    ylim([-1 1]);
end

end