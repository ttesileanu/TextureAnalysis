% make plots comparing NI predictions to measured PP thresholds for ternary
% textures

%% Load the data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Single planes

% XXX need to add the predictions
plotTernaryMatrix([], ternaryAvg, 'ellipses', false, 'mask', ~strcmp(ternaryAvg.groups, 'A_1'));