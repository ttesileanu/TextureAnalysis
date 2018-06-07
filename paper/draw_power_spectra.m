%% Choose images to focus on

images0 = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');
% images = images0(1:100);
images = images0(1);

%% Calculate power spectrum before and after whitening

filter0 = open(fullfile('filters', 'filter1x32.mat'));
filter = filter0.filter;

radial_spectra_before = cell(size(images));
radial_spectra_after = cell(size(images));
avg_fourier_before = [];
avg_fourier_after = [];
radii_before = cell(size(images));
radii_after = cell(size(images));
progress = TextProgress('power spectra');
% nbins = 24;
for i = 1:length(images)
    crt_image = loadLUMImage(fullfile('NaturalImages', images{i}));
    [crt_image_filter, crt_filter_crop] = filterImage(crt_image, filter);
    crt_image_crop = crt_image(crt_filter_crop(1):crt_filter_crop(3), crt_filter_crop(2):crt_filter_crop(4));
    
%     [radial_spectra_before{i}, radii_before{i}, crt_fourier_before] = ...
%         powerspec(crt_image_crop, 'patchSize', size(filter), 'nbins', nbins);
%     [radial_spectra_after{i}, radii_after{i}, crt_fourier_after] = ...
%         powerspec(crt_image_filter, 'patchSize', size(filter), 'nbins', nbins);
    [radial_spectra_before{i}, radii_before{i}, crt_fourier_before] = powerspec(crt_image_crop);
    [radial_spectra_after{i}, radii_after{i}, crt_fourier_after] = powerspec(crt_image_filter);
    if max(abs(radii_after{i} - radii_before{i}))
        error('The radial bins aren''t the same for before and after spectra.');
    end
    if isempty(avg_fourier_before)
        avg_fourier_before = zeros(size(crt_fourier_before));
    end
    avg_fourier_before = avg_fourier_before + crt_fourier_before;
    
    if isempty(avg_fourier_after)
        avg_fourier_after = zeros(size(crt_fourier_after));
    end
    avg_fourier_after = avg_fourier_after + crt_fourier_after;
    progress.update(100*i/length(images));
end
progress.done('done');

avg_fourier_before = avg_fourier_before / length(images);
avg_fourier_after = avg_fourier_after / length(images);

% diffs = arrayfun(@(i) max(abs(radii_before{i} - radii_before{1})), 2:length(radii_before));
%
% if max(diffs) > 1e-10
%     error('The radial bins aren''t the same for all spectra.');
% end
% 
% avg_radial_spectrum_before = mean(cell2mat(radial_spectra_before), 2);
% avg_radial_spectrum_after = mean(cell2mat(radial_spectra_after), 2);
% avg_radii_before = radii_before{1};
% avg_radii_after = radii_after{1};

[avg_radial_spectrum_before, avg_radii_before] = powerspec(avg_fourier_before, 'fromFourier', true);
[avg_radial_spectrum_after, avg_radii_after] = powerspec(avg_fourier_after, 'fromFourier', true);
if max(abs(avg_radii_before - avg_radii_after)) > 1e-10
    error('The radial bins aren''t the same for before&after spectra.');
end

%% Make plots

figSize = [1.1 0.7];

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 figSize];

coeffs = polyfit(log10(avg_radii_before(2:11)'), log10(avg_radial_spectrum_before(2:11)), 1);

loglog(avg_radii_before, 10.^(coeffs(1)*log10(avg_radii_before) + coeffs(2)), ...
    'color', [0.6, 0.6, 0.6], 'linewidth', 0.25);
hold on;
loglog(avg_radii_before, avg_radial_spectrum_before, 'color', [0    0.4470    0.7410]);

xlabel('log radius');
ylabel('log power');
set(gca, 'xtick', [], 'ytick', []);

xlim([4, 540]);
ylim([1.5e11, 4e15]);
% xlim([0.13 20]);
% ylim([5e4, 9e11]);

beautifygraph('fontscale', 0.5);
preparegraph;

safe_print(fullfile('figs', 'draft', 'power_spectrum_before_filter'));

fig = figure;
fig.Units = 'inches';
fig.Position = [3 1 figSize];

% const_y = geomean(avg_radial_spectrum_after);
const_y = 10^(trapz(log10(avg_radii_after(2:end)), log10(avg_radial_spectrum_after(2:end))) / ...
              trapz(log10(avg_radii_after(2:end)), ones(size(avg_radial_spectrum_after(2:end)))));
loglog(avg_radii_after, repmat(const_y, size(avg_radii_after)), ...
    'color', [0.6, 0.6, 0.6], 'linewidth', 0.25);
hold on;
loglog(avg_radii_after, avg_radial_spectrum_after, 'color', [0    0.4470    0.7410]);
xlabel('log radius');
ylabel('log power');
set(gca, 'xtick', [], 'ytick', []);

xlim([4, 540]);
ylim([7e7, 2e12]);

beautifygraph('fontscale', 0.5);
preparegraph;

safe_print(fullfile('figs', 'draft', 'power_spectrum_after_filter'));

%% Make 2d power spectrum plots

figSize = size(avg_fourier_before);

fig = figure;
fig.Units = 'pixels';
fig.Position = [fig.Position(1:2) fliplr(figSize)/2];

ax = axes;
ax.Position = [0 0 1 1];

data1 = log10(abs(avg_fourier_before));
% 95% of values
rng1 = [quantile(data1(:), 0.025), quantile(data1(:), 0.975)];

imagesc(data1, rng1);
colormap('gray');

axis equal;
axis off;

drawnow;
frame = getframe(fig);
imwrite(frame.cdata, fullfile('figs', 'draft', 'power_spectrum_2d_before_filter.png'));

figSize = size(avg_fourier_before);

fig = figure;
fig.Units = 'pixels';
fig.Position = [fig.Position(1:2) fliplr(figSize)/2];

ax = axes;
ax.Position = [0 0 1 1];

data2 = log10(abs(avg_fourier_after));
% keep the same dynamic range as for the 'before' spectrum
rng2 = median(data2(:)) + [-0.5 0.5]*diff(rng1);

imagesc(data2, rng2);
colormap('gray');

axis equal;
axis off;

drawnow;
frame = getframe(fig);
imwrite(frame.cdata, fullfile('figs', 'draft', 'power_spectrum_2d_after_filter.png'));