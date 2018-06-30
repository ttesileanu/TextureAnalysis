function image = imageToLUM(fname)
% imageToLUM Convert a single image to LUM Matlab format.
%   image = imageToLUM(fname) loads the image pointed to by the given
%   file name, and returns it in LUM format -- a luminance-only matrix. It
%   assumes RGB colors in sRGB profile for images with three samples per
%   pixel.

image = double(imread(fname));

if size(image, 3) == 3
    % convert to gray scale
    % this is a simple transformation that approximates the proper
    % gamma correction for sRGB
    image = sqrt(0.299*image(:, :, 1).^2 + ...
        0.587*image(:, :, 2).^2 + ...
        0.114*image(:, :, 3).^2);
elseif size(image, 3) ~= 1
    % don't really know what to do with other numbers of samples
    error([mfilename ':badsamps'], 'Can only work with 1 or 3 samples per pixel.');
end

end