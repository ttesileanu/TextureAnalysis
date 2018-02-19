% study the effects of noise on texture results

%% Generate random samples, G=3 and continuous stats

rng(3443);

G = 3;
n_samples = 1000;
patch_size = 64;
% patches = arrayfun(@(p) randi([0, G], patch_size)/(G-1), 1:n_samples, 'uniform', false);
patches = arrayfun(@(p) rand(patch_size), 1:n_samples, 'uniform', false);
evsFiniteG = [];
evsCont = [];
progress = TextProgress;
for i = 1:n_samples
    [~, crt_ev] = processBlock(quantize(patches{i}, G), G);
    [~, crt_ev_cont] = processBlock(patches{i}, inf);
    if isempty(evsFiniteG)
        evsFiniteG = zeros(n_samples, length(crt_ev));
    end
    if isempty(evsCont)
        evsCont = zeros(n_samples, length(crt_ev_cont));
    end
    evsFiniteG(i, :) = crt_ev; %#ok<SAGROW>
    evsCont(i, :) = crt_ev_cont; %#ok<SAGROW>
    if mod(i, 10) == 0
        progress.update(100*i/n_samples);
    end
end
progress.done;

%% Shape of noise around origin, G=3 stats

% keep things reproducible
rng(32423);

% find all group names
mtc = processBlock('mtc', G);
groups = mtc.coord_groups;
directions = cell(1, G-1);
for i = 1:G-1
    crt_dir = num2str((1:G) == i, '%d');
    directions{i} = ['[' crt_dir ']'];
end

% generate labels for all the coordinates in ev
labels = cell(1, length(groups)*(G-1));
for i = 1:length(labels)
    crt_group_idx = floor((i-1) / (G-1)) + 1;
    crt_dir_idx = mod(i-1, G-1) + 1;
    labels{i} = [groups{crt_group_idx}.name directions{crt_dir_idx}];
end

% choose how many plots to make
n_coords = size(evsFiniteG, 2);
n_plots = n_coords;

% set a range for the histograms
group_stds = arrayfun(@(i) std(evsFiniteG(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsFiniteG(:, i)) + hist_bins;
    histogram(evsFiniteG(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
end

%% Shape of noise around origin, continuous stats

% keep things reproducible
rng(32423);

% find all group names
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% choose how many plots to make
n_coords = size(evsCont, 2);
n_plots = n_coords;

% choose a range for the histograms
group_stds = arrayfun(@(i) std(evsCont(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsCont(:, i)) + hist_bins;
    histogram(evsCont(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
end

%% Susceptibility to noise away from origin

rng(23498);

patch_size = 64;
generator = PatchAxisGenerator('AB_1_1', [0 0 1], patch_size);
% generator = PatchAxisGenerator('AB_1_2', [0 0 1], patch_size);
generator.locations = 0.5;
generator.next;
base_patch = generator.samples;

G = 3;
n_samples = 1000;
noise_size = 0.3;
patchesFar = arrayfun(@(p) min(max(base_patch + noise_size*(2*rand(patch_size)-1), 0), 1), ...
    1:n_samples, 'uniform', false);
evsFarFiniteG = [];
evsFarCont = [];
progress = TextProgress;
for i = 1:n_samples
    [~, crt_ev] = processBlock(quantize(patchesFar{i}, G), G);
    [~, crt_ev_cont] = processBlock(patchesFar{i}, inf);
    if isempty(evsFarFiniteG)
        evsFarFiniteG = zeros(n_samples, length(crt_ev));
    end
    if isempty(evsFarCont)
        evsFarCont = zeros(n_samples, length(crt_ev_cont));
    end
    evsFarFiniteG(i, :) = crt_ev; %#ok<SAGROW>
    evsFarCont(i, :) = crt_ev_cont; %#ok<SAGROW>
    if mod(i, 10) == 0
        progress.update(100*i/n_samples);
    end
end
progress.done;

%% Shape of noise far from origin, G=3 stats

% keep things reproducible
rng(32423);

% find all group names
mtc = processBlock('mtc', G);
groups = mtc.coord_groups;
directions = cell(1, G-1);
for i = 1:G-1
    crt_dir = num2str((1:G) == i, '%d');
    directions{i} = ['[' crt_dir ']'];
end

% generate labels for all the coordinates in ev
labels = cell(1, length(groups)*(G-1));
for i = 1:length(labels)
    crt_group_idx = floor((i-1) / (G-1)) + 1;
    crt_dir_idx = mod(i-1, G-1) + 1;
    labels{i} = [groups{crt_group_idx}.name directions{crt_dir_idx}];
end

% choose how many plots to make
n_coords = size(evsFarFiniteG, 2);
n_plots = n_coords;

% set a range for the histograms
group_stds = arrayfun(@(i) std(evsFarFiniteG(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsFarFiniteG(:, i)) + hist_bins;
    histogram(evsFarFiniteG(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
end

%% Shape of noise far from origin, continuous stats

% keep things reproducible
rng(32423);

% find all group names
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% choose how many plots to make
n_coords = size(evsFarCont, 2);
n_plots = n_coords;

% set a range for the histograms
group_stds = arrayfun(@(i) std(evsFarCont(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsFarCont(:, i)) + hist_bins;
    histogram(evsFarCont(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
end