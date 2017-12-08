 function showResultsStats(res, varargin)
% showResultsStats Display distributions of texture statistics.
%   showResultsStats(res) makes scatter plots showing the distribution of
%   texture statistics for a set of results.
%
%   Additional arguments are directly passed to smartscatter.
%
%   See also: SMARTSCATTER.

%figure;
set(gcf, 'position', [50 50 961 600], 'paperposition', [0.25 0.25 9.61 6]);

% get rid of spurious patches
ev = res.ev;
mask = (ev(:, 1) ~= 1);

c_labels = {'\gamma', '\beta_|', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{|-}', '\theta_{-|}', '\theta_{\_|}', '\theta_{|\_}', '\alpha'};

order = [1, 10, 2, 3, 4, 5, 6, 7, 8, 9];

allObjIds = unique(res.objIds);
colors = parula(max(allObjIds)+1);

for i0 = 1:2:10
    i1 = order(i0);
    i2 = order(i0+1);
    
    subplot(2, 3, (i0+1)/2);
    hold on;
    for objId = 1:max(allObjIds)
        crtMask = (mask & res.objIds == objId);
        color = colors(objId+1, :);
        if sum(crtMask) > 0
            smartscatter(ev(crtMask, i1), ev(crtMask, i2), 10, color, 'filled', varargin{:});
        end
    end
    
    xlabel(['comp. ' int2str(i1) ': ' c_labels{i1}]);
    ylabel(['comp. ' int2str(i2) ': ' c_labels{i2}]);    
end

end