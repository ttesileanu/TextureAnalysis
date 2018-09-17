function m = ternaryPCToRGB(pcs)
% ternaryPCToRGB Convert principal component vectors for ternary textures
% to RGB colors.
%   m = ternaryPCToRGB(pcs) converts the principal components for ternary
%   textures that are organized along the columns of the `pcs` argument to
%   a 3d array of RGB colors, such that `m(i, j, :)` represents the color
%   for texture plane `j` in the principal component `i`. The input is
%   assumed to be in the full-probability, 99-dimensional representation in
%   which each consecutive group of 3 rows of `pcs` corresponds to one
%   texture plane.

% check input
if mod(size(pcs, 1), 3) ~= 0
    error([mfilename ':badsz'], 'The input should have a number of rows that is a multiple of 3.');
end

% allocate output structure
nPlanes = size(pcs, 1) / 3;
nVecs = size(pcs, 2);

% m = zeros(nVecs, nPlanes, 3);

% start by converting 3d probabilities to luminances (departure from 0) and
% normalized 3-vectors
all3d = reshape(pcs, 3, []);
luminance = sqrt(sum(all3d.^2, 1));
% set things below a certain tolerance to 0
tol = 1e-4;
luminance(abs(luminance) < tol) = 0;
% normalize so that each non-zero 3-component vector has L2 norm 1
mask = (luminance > 0);
all3dnorm(:, mask) = bsxfun(@rdivide, all3d(:, mask), luminance(mask));
% ...and clip other vectors to 0
all3dnorm(:, ~mask) = 0;

% make this a 3-column matrix
all3dnorm = all3dnorm';

% saturation = std(all3dnorm, [], 1);

% find entry that has a sign different from the others
signs = sign(all3dnorm);
% count 0 as "-" so that when we have a 0, a +, and a -, we choose the +
signs(signs == 0) = -1;
specialSigns = (signs ~= sum(signs, 2));

% values that determine colors are taken from the special signs
colorDeterminants = all3dnorm(specialSigns);

colorBasis = {[0.7 0 0 ; 0 0 1], [0.55 0.55 0 ; 0 0.55 0.55], [0 0.5 0 ; 0.45 0 0.45]};
colors = zeros(size(all3dnorm));
for i = 1:3
    mask = specialSigns(:, i);
    crtBasis = colorBasis{i};
    colors(mask, :) = 0.5*colorDeterminants(mask)*(crtBasis(1, :) - crtBasis(2, :)) + ...
        0.5*(crtBasis(1, :) + crtBasis(2, :));
end

% multiply the luminance back in
colors = bsxfun(@times, colors, luminance(:));

% hue = (2 + all3dnorm(1, :) - all3dnorm(3, :))/4;
% 
% % normalize luminances so that maximum is 1 in every principal component
% luminanceMatrix = reshape(luminance, [nPlanes, nVecs]);
% luminanceMatrix = bsxfun(@rdivide, luminanceMatrix, max(luminanceMatrix, [], 1));
% luminance = luminanceMatrix(:);
% 
% % generate RGB colors
% rgb = hsv2rgb(max(min([hue(:) saturation(:) luminance(:)], 1), 0));

m = permute(reshape(colors, [nPlanes, nVecs, 3]), [2 1 3]);

end
