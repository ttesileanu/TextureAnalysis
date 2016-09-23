function image = convertToLog(image, threshold)
% convertToLog Convert image to logarithmic space.
%   newImage = convertToLog(image) converts the input matrix to logarithmic
%   space by taking the log of all elements. Elements smaller than 1e-6 are
%   clipped at that value.
%
%   convertToLog(image, threshold) uses the given threshold instead of 1e-6.

if nargin < 2
    threshold = 1e-6;
end

image(image < threshold) = threshold;
image = log(image);

end