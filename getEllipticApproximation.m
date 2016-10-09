function [center, sigmas] = getEllipticApproximation(mask)
% getEllipticApproximation Find center and horizontal and vertical standard
% deviations for a region in an image.
%   [center, sigmas] = getEllipticApproximation(mask) calculates the
%   centroid of the 'on' pixels in the logical image `mask`. It also
%   calculates the extent of the mask on the x- and y-axes.

n = sum(mask(:));

if n == 0
    center = nan(1, 2);
    sigmas = nan(1, 2);
    return;
end

center = zeros(1, 2);
[yGrid, xGrid] = ndgrid(1:size(mask, 1), 1:size(mask, 2));
center(1) = sum(sum(xGrid.*mask))/n;
center(2) = sum(sum(yGrid.*mask))/n;

sigmas = zeros(1, 2);
sigmas(1) = sqrt(sum(sum(xGrid.^2 .* mask))/n - center(1)^2);
sigmas(2) = sqrt(sum(sum(yGrid.^2 .* mask))/n - center(2)^2);

end