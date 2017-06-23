function sharp = getSharpnessMeasure(image)
% getSharpnessMeasure Get a measure of the sharpness of a grayscale image.
%   sharp = getSharpnessMeasure(image) returns a measure of the sharpness
%   of the input image. This is done by convolving with a Laplacian filter
%   and calculating the exponential of the median of the logarithm of the
%   absolute value of the result, after the image is normalized so that the
%   median pixel value is 1. NaN pixels in the filtered image are ignored.

medianValue = median(image(:), 'omitnan');
if medianValue <= 0
    sharp = 0;
else 
    image = image/medianValue;
    lapFlt = [0 1 0; 1 -4 1; 0 1 0];
    filtered = conv2(image, lapFlt, 'valid');
%    filtered = imgradient(image);
%    sharp = median(abs(filtered(:)), 'omitnan');
    sharp = exp(median(log(abs(filtered(:))), 'omitnan'));
end

end