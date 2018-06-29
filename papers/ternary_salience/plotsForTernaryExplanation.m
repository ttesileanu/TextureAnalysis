% make plots that will be used to explain the ternary texture analysis

%% Save independent texture patches

% this will be used on the tips and in the center of a ternary alloy plot

% keep things reproducible
rng(3434);

% choose the plane, locations in each plane (with labels for saving), and patch size
planes = {'AB_1_2'};
texAxes = {[0 1 0], [0 0 1], [1 0 0], [1/3, 1/3, 1/3]};
locLabels = {'1', '2', '0', 'c'};
patchSize = 128;

% draw everything
for i = 1:length(planes)
    plane = planes{i};
    
    for k = 1:length(texAxes)
        crtLoc = texAxes{k};
        
        % generate a texture all the way at the end of the texture axis
        generator = PatchAxisGenerator(plane, texAxes{k}, patchSize);
        generator.locations = 1;
        generator.next;
        crtPatch = generator.samples;
        
        % draw the generated patch
        fig = figure;
        fig.Units = 'pixels';
        fig.Position = [fig.Position(1:2) patchSize patchSize];
        
        ax = axes;
        ax.Position = [0 0 1 1];
        imagesc(crtPatch, [0 1]);
        colormap('gray');
        
        axis equal;
        axis off;
        
        % save the frame data as a PNG
        drawnow;
        frame = getframe(fig);
        imwrite(frame.cdata, fullfile('figs', 'draft', ...
            ['ternaryPatch_' plane '_' locLabels{k} '.png']));
    end
end

%% Draw simplex and guiding circles

% where to draw inset bar plots
texAxes = {[0 1 0], [0 0 1], [1 0 0]};
% color of bars
barColor = [64, 90, 128]/255;
% size of bar plot
patchDim = 0.5;
% shift of bar plot compared to axis coordinates
patchShift = 0.07;

% build the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2 2];

hold on;
axis equal;

% draw the guide simplex and circles
drawTernaryTriangle('edgelabels', 'none');

% draw
for k = 1:length(texAxes)
    crtLoc = texAxes{k};
    
    % project 3d axes coordinates to 2d alloy plot coordinates
    crtT = ternary3to2(crtLoc);
    
    if k == 3
        crtT(1) = crtT(1) - 4*patchShift;
    elseif k == 2
        crtT(2) = crtT(2) - patchDim;
        crtT(1) = crtT(1) - 4*patchShift;
    else
        crtT(2) = crtT(2) - patchDim;
    end
    
    % draw the bar plots
    % the weird ifs just position the plots in more-or-less reasonable places
    if k ~= 3
        insetbar([crtT + [0 1.15*patchDim] patchDim patchDim/2], crtLoc + 0.05, ...
            'color', barColor);
    else
        insetbar([crtT + [0 -0.65*patchDim] patchDim patchDim/2], crtLoc + 0.05, ...
            'color', barColor);
    end
end

% draw a bar plot in the center
insetbar([0 patchShift patchDim patchDim/2], [1/3 1/3 1/3] + 0.05, ...
    'color', barColor);

% update the looks of the axes
beautifygraph('ticks', 'off', 'ticklabels', false, ...
    'titlesize', 10, 'titleweight', 'normal', 'noaxes', true, ...
    'fontscale', 0.9);
preparegraph;

safePrint(fullfile('figs', 'draft', 'exampleTernaryPlane'));