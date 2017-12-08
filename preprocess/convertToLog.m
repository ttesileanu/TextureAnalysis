function image = convertToLog(image, threshold, type)
% convertToLog Convert image to logarithmic space.
%   newImage = convertToLog(image) converts the input matrix to logarithmic
%   space by taking the log of all elements. Elements smaller than 1e-6 are
%   clipped at that value.
%
%   convertToLog(image, threshold) uses the given threshold instead of 1e-6.

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

end