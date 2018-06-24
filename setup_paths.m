function setup_paths
% SETUP_PATHS Temporarily add texture analysis folders to Matlab's path.
%   SETUP_PATHS adds all the folders containing texture analysis routines
%   to Matlab's path. The change will only stay in effect until Matlab
%   quits, so SETUP_PATHS needs to be called every time Matlab is
%   restarted.

% add current folder
PWD = pwd;
addpath(PWD);

% add subfolders
addpath(fullfile(PWD, 'analysis'));
addpath(fullfile(PWD, 'batch'));
addpath(fullfile(PWD, 'filters'));
addpath(fullfile(PWD, 'graphing'));
%addpath(fullfile(PWD, 'gui'));
addpath(fullfile(PWD, 'indices'));
addpath(fullfile(PWD, 'load'));
addpath(fullfile(PWD, 'misc'));
addpath(fullfile(PWD, 'patch'));
addpath(fullfile(PWD, 'preprocess'));
addpath(fullfile(PWD, 'postprocess'));
%addpath(fullfile(PWD, 'sandbox'));
%addpath(fullfile(PWD, 'segment'));

end