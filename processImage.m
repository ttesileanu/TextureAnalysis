function img = processImage(fname)
% processImage Convert a single image to LUM Matlab format.
%   image = processImage(fname) loads the image pointed to by the given
%   file name, and returns it in LUM format -- a luminance-only matrix.

img = double(imread(fname));

if size(img, 3) == 3
    % convert to gray scale
    % this is a simple transformation that approximates the proper
    % gamma correction for sRGB
    img = sqrt(0.299*img(:, :, 1).^2 + ...
        0.587*img(:, :, 2).^2 + ...
        0.114*img(:, :, 3).^2);
elseif size(img, 3) ~= 1
    % don't really know what to do with other numbers of samples
    error([mfilename ':badsamps'], 'Can only work with 1 or 3 samples per pixel.');
end

end