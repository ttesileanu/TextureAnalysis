function [image, crop] = logTransform(image, threshold, type)
% logTransform Convert image to logarithmic space.
%   newImage = logTransform(image) converts the input matrix to logarithmic
%   space by taking the log of all elements. Elements smaller than 1e-6 are
%   clipped at that value before taking the log.
%
%   logTransform(image, threshold) uses the given threshold instead of 1e-6.
%
%   [newImage, crop] = logTransform(...) returns a crop for the log
%   transform. This is always equal to [1 1 size(image)].

if nargin < 2
    threshold = 1e-6;
end
if nargin < 3
    type = '';
end

% image(image < threshold) = threshold;
if isempty(type)
    image = max(image, threshold);
elseif strcmp(type, 'legacy')
    image(image <= 0) = threshold;
else
    error([mfilename ':badtype'], 'Unknown type.');
end

image = log(image);
crop = [1 1 size(image)];

end