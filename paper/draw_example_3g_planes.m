% make diagrams of some sample 3-gray-level planes

%% setup

% keep this reproducible
rng(3434);

planes = {'AB_1_2'};
tex_axes = {[0 1 0], [0 0 1], [1 0 0]};
patch_size = 16;
patch_dim = 0.5;
patch_shift = 0.07;
bar_color = [64, 90, 128]/255;

for i = 1:length(planes)
    plane = planes{i};
    
    plane_order = find(plane == '_', 1) - 1;
    plane_eq = '';
    for k = 1:plane_order
        crt_coeff = str2double(plane(plane_order + 2*k));
        if crt_coeff > 1
            plane_eq = [plane_eq int2str(crt_coeff)]; %#ok<AGROW>
        end
        plane_eq = [plane_eq plane(k)]; %#ok<AGROW>
        if k < plane_order
            plane_eq = [plane_eq '+']; %#ok<AGROW>
        end
    end
    
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [fig.Position(1:2) 2 2];
    
    hold on;
    
    drawTernaryTriangle('edgelabels', 'none');
    
    max_t2 = sqrt(3)/2;
%     angle_range = linspace(0, 2*pi, 100);
%     
%     % draw circles for orientation (radius 1 and 1/2)
%     plot(cos(angle_range), sin(angle_range), ':', 'color', [0.4 0.4 0.4]);
%     plot(0.5*cos(angle_range), 0.5*sin(angle_range), ':', 'color', [0.4 0.4 0.4]);
%     
%     % draw the probability triangle
%     plot([-1/2 1 -1/2 -1/2], [-max_t2 0 max_t2 -max_t2], 'color', [0.5 0.7 1]);
%     
%     % draw the main axes
%     plot([0 1.5], [0 0], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
%     plot([0 -1/2*1.5], [0 1.5*max_t2], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
%     plot([0 -1/2*1.5], [0 -1.5*max_t2], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
%     
%     % label the corners
%     text(0.95, -0.1, '[0,1,0]', 'fontsize', 12);
%     text(0.82, -0.25, [plane_eq '=1 (mod 3)'], 'fontsize', 12);
%     
%     text(-0.95,  max_t2+0.05, '[0,0,1]', 'fontsize', 12);
%     text(-1.35,  max_t2-0.10, [plane_eq '=2 (mod 3)'], 'fontsize', 12);
%     
%     text(-0.95, -max_t2-0.01, '[1,0,0]', 'fontsize', 12);
%     text(-1.35, -max_t2-0.16, [plane_eq '=0 (mod 3)'], 'fontsize', 12);
    
    axis equal;
    
    for k = 1:length(tex_axes)
        crt_loc = tex_axes{k};
        
        % this is a litle roundabout, but easier conceptually
        crt_t1 = (3*crt_loc(2) - 1)/2;
        crt_t2 = (crt_loc(3) - crt_loc(1)) * max_t2;
        
        generator = PatchAxisGenerator(plane, tex_axes{k}, patch_size);
        generator.locations = 1;
        generator.next;
        crt_patch = generator.samples;
        
        colormap('gray');
        n_cmap = size(colormap, 1);
        edge = patch_dim/patch_size/2;

        if k == 3
            crt_t2 = crt_t2 - patch_dim - patch_shift;
            crt_t1 = crt_t1 - patch_dim + patch_shift;
        elseif k == 2
            crt_t2 = crt_t2 + patch_shift;
            crt_t1 = crt_t1 - patch_dim + patch_shift;
        else
            crt_t2 = crt_t2 + patch_shift;
            crt_t1 = crt_t1 + patch_shift;
        end
        
        image([crt_t1 crt_t1 + patch_dim], [crt_t2 crt_t2 + patch_dim], floor(crt_patch*n_cmap));
        plot([crt_t1-edge crt_t1+patch_dim+edge crt_t1+patch_dim+edge crt_t1-edge crt_t1-edge], ...
             [crt_t2-edge crt_t2-edge crt_t2+patch_dim+edge crt_t2+patch_dim+edge crt_t2-edge], ...
             'k');
        
        if k ~= 3
            insetbar([crt_t1 crt_t2+1.15*patch_dim patch_dim patch_dim/2], crt_loc + 0.05, ...
                'color', bar_color);
        else
            insetbar([crt_t1 crt_t2-0.65*patch_dim patch_dim patch_dim/2], crt_loc + 0.05, ...
                'color', bar_color);
        end
    end
    
    % draw a patch in the center
    crt_t1 = patch_shift;
    crt_t2 = patch_shift;
    image([crt_t1 crt_t1+patch_dim], [crt_t2 crt_t2+patch_dim], n_cmap/2*randi([0 2], patch_size));
    plot([crt_t1-edge crt_t1+patch_dim+edge crt_t1+patch_dim+edge crt_t1-edge crt_t1-edge], ...
             [crt_t2-edge crt_t2-edge crt_t2+patch_dim+edge crt_t2+patch_dim+edge crt_t2-edge], ...
             'k');
    insetbar([crt_t1 crt_t2+1.15*patch_dim patch_dim patch_dim/2], [1/3 1/3 1/3] + 0.05, ...
        'color', bar_color);
    
%    h = title(plane);
%    h.Position = [h.Position(1) h.Position(2) - 0.1];
    
    beautifygraph('ticks', 'off', 'ticklabels', false, ...
        'titlesize', 10, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.9);
    preparegraph;
    
    safe_print(fullfile('figs', 'draft', ['example_3g_plane_' plane]));
end

%% Save independent patches

rng(3434);

planes = {'AB_1_2'};
tex_axes = {[0 1 0], [0 0 1], [1 0 0], [1/3, 1/3, 1/3]};
loc_labels = {'1', '2', '0', 'c'};
patch_size = 128;

for i = 1:length(planes)
    plane = planes{i};
    
    for k = 1:length(tex_axes)
        crt_loc = tex_axes{k};
        
        generator = PatchAxisGenerator(plane, tex_axes{k}, patch_size);
        generator.locations = 1;
        generator.next;
        crt_patch = generator.samples;
        
        fig = figure;
        fig.Units = 'pixels';
        fig.Position = [fig.Position(1:2) patch_size patch_size];
        
        ax = axes;
        ax.Position = [0 0 1 1];
        imagesc(crt_patch, [0 1]);
        colormap('gray');
        
        axis equal;
        axis off;
        
        drawnow;
        frame = getframe(fig);
        imwrite(frame.cdata, fullfile('figs', 'draft', ...
            ['patch_' plane '_' loc_labels{k} '.png']));
    end
end

%% Draw simplex and guiding circles

% keep this reproducible
rng(3434);

tex_axes = {[0 1 0], [0 0 1], [1 0 0]};
bar_color = [64, 90, 128]/255;
patch_dim = 0.5;
patch_shift = 0.07;

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2 2];

hold on;

drawTernaryTriangle('edgelabels', 'none');

max_t2 = sqrt(3)/2;

axis equal;

for k = 1:length(tex_axes)
    crt_loc = tex_axes{k};
    
    % this is a litle roundabout, but easier conceptually
    crt_t1 = (3*crt_loc(2) - 1)/2;
    crt_t2 = (crt_loc(3) - crt_loc(1)) * max_t2;
    
    if k == 3
        crt_t1 = crt_t1 - 4*patch_shift;
    elseif k == 2
        crt_t2 = crt_t2 - patch_dim;
        crt_t1 = crt_t1 - 4*patch_shift;
    else
        crt_t2 = crt_t2 - patch_dim;
        %             crt_t1 = crt_t1 + patch_shift;
    end
    
    if k ~= 3
        insetbar([crt_t1 crt_t2+1.15*patch_dim patch_dim patch_dim/2], crt_loc + 0.05, ...
            'color', bar_color);
    else
        insetbar([crt_t1 crt_t2-0.65*patch_dim patch_dim patch_dim/2], crt_loc + 0.05, ...
            'color', bar_color);
    end
end

% draw a patch in the center
%     crt_t1 = patch_shift;
%     crt_t2 = patch_shift;
crt_t1 = 0;
crt_t2 = 0;
insetbar([crt_t1 crt_t2 + patch_shift patch_dim patch_dim/2], [1/3 1/3 1/3] + 0.05, ...
    'color', bar_color);

beautifygraph('ticks', 'off', 'ticklabels', false, ...
    'titlesize', 10, 'titleweight', 'normal', 'noaxes', true, ...
    'fontscale', 0.9);
preparegraph;

safe_print(fullfile('figs', 'draft', 'example_3g_plane_background'));