function scatterTextureStats(ev, varargin)
% scatterTextureStats Make scatter plots of texture statistics.
%   scatterTextureStats(ev) makes a figure containing 5 scatter plots
%   showing the texture statistics in the input array `ev`.
%
%   Additional arguments are directly passed to smartscatter.
%
%   scatterTextureStats({ev1, ev2, ...}) makes a plot containing several
%   sets of texture statistics. They are colored using the parula color
%   map.
%
%   scatterTextureStats({{ev1, ...}, [c1, ...]) uses a coloring scheme
%   based on the values c1, ... Specifically, c == 0 is mapped to black
%   (the middle of the color map), c == 1 to red (the end), and c == -1 to
%   blue (the beginning).
%
%   See also: smartscatter.

labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

fig = figure;
fig.Units = 'inches';
oldPos = fig.Position;
centroid = oldPos(1:2) + 0.5*oldPos(3:4);
figSize = [10 1.6];
fig.Position = [centroid - 0.5*figSize figSize];

pairs = {[1 10], [2 3], [4 5], [6, 7], [8 9]};

if iscell(ev{1}) && numel(ev) == 2 && isnumeric(ev{2})
    has_color_scheme = true;
else
    has_color_scheme = false;
end

ax = cell(1, length(pairs));
for i = 1:length(pairs)
    idx1 = pairs{i}(1);
    idx2 = pairs{i}(2);
    
    ax{i} = axes;
    ax{i}.Units = 'normalized';
    ax{i}.OuterPosition = [(i-1)/length(pairs) 0 1/length(pairs) 0.9];
    
%    subplot(1, length(pairs), i);
    if ~iscell(ev)
        smartscatter(ev(:, idx1), ev(:, idx2), varargin{:});
    else
        if has_color_scheme
            colors0 = redbluem(255);
            cidxs = 1 + round(0.5*(1 + ev{2})*(size(colors0, 1) - 1));
            colors = colors0(cidxs, :);
            ev_here = ev{1};
        else
            colors = redbluem(length(ev));
            ev_here = ev;
        end
        hold on;
        for j = 1:length(ev_here)
            crt_ev = ev_here{j};
            smartscatter(crt_ev(:, idx1), crt_ev(:, idx2), [], colors(j, :), varargin{:});
        end
    end
    xlabel(labels{idx1});
    ylabel(labels{idx2});
    beautifygraph;
end
mean_bottom = mean(cellfun(@(a) a.Position(2), ax));
mean_height = mean(cellfun(@(a) a.Position(4), ax));
for i = 1:length(ax)
    ax{i}.Position(2) = mean_bottom;
    ax{i}.Position(4) = mean_height;
end

preparegraph;

end

function cmap = redbluem(n)
% Generate a red-blue colormap that goes through black instead of white.

cmap = redblue(n);
cmap = max(min(cmap - (sum(cmap, 2) - 1)/2, 1), 0);

end