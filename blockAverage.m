function out = blockAverage(img, n, type)
% blockAverage Scale down an image using block averaging.
%   out = blockAverage(img, n) returns an image that was scaled down
%   by a factor `n` using block averaging.
%
%   out = blockAverage(img, n, type) uses an averaging method given by
%   `type`. Allowed values for `type` are:
%       'avg' (default)
%           The image is down-sized by averaging over n x n blocks.
%       'sub'
%           Naive subsampling: every nth pixel is kept.
%       'wta'
%           This is 'winner takes all': like 'avg', but returns a logical
%           matrix which is false whenever the block-averaged output is
%           smaller than 0.5, and true otherwise.

if nargin < 3
    type = 'avg';
end

% check the input
if ~ismember(type, {'avg', 'wta', 'sub'})
    error([mfilename ':badtype'], 'Unrecognized averaging type.');
end

% to average over n x n blocks, convolve with constant matrix before
% subsampling
if ismember(type, {'avg', 'wta'})
    kernel = ones(n) / n^2;
    img = conv2(img, kernel, 'valid');
end

% we always subsample
out = img(1:n:end, 1:n:end);

% for winner takes all, just threshold
if strcmp(type, 'wta')
    out = (out >= 0.5);
end

end