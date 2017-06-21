% perform a multiscale analysis

%% Generate multiscale filters

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

originalR = 60;
N_values = 1:6;

filters = cell(1, length(N_values)); %#ok<*NASGU>
for i = 1:length(N_values)
    crtN = N_values(i);
    filters{i} = generateFourierWhitenFilter(images, 'NaturalImages', crtN, round(originalR/crtN));
end

%% Load multiscale filters

originalR = 60;
N_values = 1:6;

filters = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtFilter = open(['filter' int2str(crtN) '.mat']);
    filters{i} = crtFilter.ans;
end

%% Generate the multiscale stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    res{i} = analyzeImageSet(images, 'NaturalImages', crtN, filters{i}, inf);
end

%% Load the multiscale stats

load('natural_nosky_multiscale.mat');
res = {res1 res2 res3 res4 res5 res6};
clear res1 res2 res3 res4 res5 res6;

%% Run cross-scale PCA

% either including all scales, or only AF=2 and up
all_ev = [];
all_ev_ge2 = [];
for i = 1:length(res)
    all_ev = [all_ev res{i}.ev]; %#ok<AGROW>
    if N_values(i) > 1
        all_ev_ge2 = [all_ev_ge2 res{i}.ev]; %#ok<AGROW>
    end
end

[pcs, pc_proj, pc_vars, ~, pc_explained] = pca(all_ev);
[pcs_ge2, pc_proj_ge2, pc_vars_ge2, ~, pc_explained_ge2] = pca(all_ev_ge2);

%% Run per-scale PCA

each_pc = cell(1, length(res));
each_pc_proj = cell(1, length(res));
each_pc_var = cell(1, length(res));
each_pc_explained = cell(1, length(res));
for i = 1:length(res)
    [each_pc{i}, each_pc_proj{i}, each_pc_var{i}, ~, each_pc_explained{i}] = ...
        pca(res{i}.ev);
end

%% Make Scree plots

fig = figure;
fig.Units = 'inches';
fig.Position = [4 0.5 8 8];

% draw cross-scale plot, AF1--6
ax_cross = axes;
ax_cross.Units = 'normalized';
ax_cross.OuterPosition = [0 0.6 0.5 0.4];

plot(pc_explained, '.-');
hold on;
plot(cumsum(pc_explained), '.-');
grid on;

title('Cross-scale PCA, N=1 to N=6, N*R=60');
legend({'per PC', 'cumulative'});

xlabel('PC#');
ylabel('Variance explained');

ylim([0 100]);

beautifygraph;

% draw cross-scale plot excluding AF1
ax_cross_ge2 = axes;
ax_cross_ge2.Units = 'normalized';
ax_cross_ge2.OuterPosition = [0.5 0.6 0.5 0.4];

plot(pc_explained_ge2, '.-');
hold on;
plot(cumsum(pc_explained_ge2), '.-');
grid on;

title('Cross-scale PCA, N=2 to N=6, N*R=60');
legend({'per PC', 'cumulative'});

xlabel('PC#');
ylabel('Variance explained');

ylim([0 100]);

beautifygraph;

% draw within-scale PCA plots
npx = ceil(sqrt(length(res)));
npy = ceil(length(res) / npx);
ax_separated = cell(length(res), 1);
for i = 1:length(res)
    crt_x = mod(i-1, npx);
    crt_y = floor((i-1) / npx);
    
    ax_separated{i} = axes;
    ax_separated{i}.Units = 'normalized';
    % crt_y == 0     --> ax_cross.OuterPosition(2)*(1 - 1/npy)
    % crt_y == npy-1 --> 0
    ax_separated{i}.OuterPosition = [...
        crt_x/npx ax_cross.OuterPosition(2)*(1 - 1/npy)*(1 - crt_y/(npy-1)) ...
        1/npx ax_cross.OuterPosition(2)/npy];

%    subplot(npy, npx, i);
    plot(each_pc_explained{i}, '.-');
    
    hold on;
    plot(cumsum(each_pc_explained{i}), '.-');
    grid on;
    
    title(['Scree plot N=' int2str(N_values(i)) ', R=' int2str(round(originalR/N_values(i)))]);
    legend({'per PC', 'cumulative'}, 'location', 'east');
    
    xlabel('PC#');
    ylabel('Variance explained');
    
    ylim([0 100]);
    
    beautifygraph;
end

preparegraph;

%% Make PC plot, cross-scale AF1 to AF6 

makeMultiscalePCAPlot(pcs, pc_explained, N_values);

%% Make PC plot,  cross-scale AF2 to AF6 

makeMultiscalePCAPlot(pcs_ge2, pc_explained_ge2, N_values(2:end));

%% Make PC plot, per scale

makeMultiscalePCAPlot(each_pc, each_pc_explained, N_values);