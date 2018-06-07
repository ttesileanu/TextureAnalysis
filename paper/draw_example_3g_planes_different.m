% make diagrams of some sample 3-gray-level planes

%% setup

% keep this reproducible
rng(3434);

planes = {'AB_1_2'};
% generate small patches on a grid
patch_size = 6;
patch_dim = 0.1;

for i = 1:length(planes)
    % prepare the figure for this plane
    plane = planes{i};
    
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [fig.Position(1:2) 2 2];
    
    hold on;
    
    axis equal;
    
    colormap('gray');
    n_cmap = size(colormap, 1);
    
    % sweep through the grid and generate patches
    for x = -0.5:patch_dim:1
        for y = -1:patch_dim:1
            v = ternary2to3([x, y]);
            % skip locations that are outside the probability simplex
            if any(v) < 0
                continue;
            end
            
            % generate the patch
            generator = PatchAxisGenerator(plane, v, patch_size);
            generator.locations = 1;
            generator.next;
            crt_patch = generator.samples;
            
            % draw the patch
            image([x-patch_dim/2, x+patch_dim/2], [y-patch_dim/2, y+patch_dim/2], ...
                floor(crt_patch*n_cmap));
                
        end
    end
    
    drawTernaryTriangle('edgelabels', 'none');
    
%    h = title(plane);
%    h.Position = [h.Position(1) h.Position(2) - 0.1];
    
    beautifygraph('ticks', 'off', 'ticklabels', false, ...
        'titlesize', 10, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.9);
    preparegraph;
    
%     safe_print(fullfile('figs', 'draft', ['example_3g_plane_' plane]));
end