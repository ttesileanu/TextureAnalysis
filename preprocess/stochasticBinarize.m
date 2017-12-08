function img = stochasticBinarize(img)
% stochasticBinarize Binarize images using a stochastic approach.
%   imgOut = stochasticBinarize(img) binarizes the matrix `img` by using a
%   stochastic approach: a pixel will be turned white (i.e., equal to 1)
%   with probability given by the initial value of the pixel, `img(i, j)`.

img = 1*(rand(size(img)) < img);

end