function [out] = blockAverage(img,n,type)
% blockAverage Scale down an image using block averaging.
%   out = blockAverage(img, n) returns an image that was scaled down
%   by a factor 'n' using block averaging.
%
%   out = blockAverage(img, n, type) uses an averaging method given by
%   'type'. Note that the factor by which the image is scaled down can be
%   2^n in this case, instead of n (when 'type' is 'sub'). Allowed values
%   for 'type' are:
%       'avg' (default)
%           The image is down-sized by averaging over n x n blocks.
%       'sub'
%           Naive subsampling by a factor of 2^(n-1). Every (2^(n-1))th
%           pixel is kept.
%       'wta'
%           This is 'winner takes all': like 'avg', but returns a logical
%           matrix which is false whenever the block-averaged output is
%           smaller than 0.5, and true otherwise.

if nargin<3
    type = 'avg';
end

if (strcmpi(type,'sub'))
    out = img(1:2^(n-1):end,1:2^(n-1):end);
    return;
end

C=ones(n)/(n^2);
outN=floor(size(img)/n);
R=conv2(img,C,'valid');
out=R( 1 + n*(0:(outN(1)-1)), 1 + n*(0:(outN(2)-1)) );

if strcmp(type,'wta')
    out = out >= 0.5;
end

end