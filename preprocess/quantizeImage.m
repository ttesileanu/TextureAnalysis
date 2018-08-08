function [img, crop] = quantizeImage(img, n, varargin)
% quantizeImage Quantize images to discrete color levels.
%   imgOut = quantizeImage(img, n) discretizes the matrix `img` to `n`
%   levels, equally spaced between 0 and 1. Values outside the [0, 1] range
%   are clipped to 0 or 1.
%
%   [imgOut, crop] = quantizeImage(...) returns a cropping area, always
%   equal to [1 1 size(img)] if `img` is a matrix, and [] otherwise.
%
%   Options:
%    'cutoffs'
%       When provided, this should be an (n-1)-element vector giving the
%       cut-off points used for discretization. If not provided, equally
%       spaced levels are used. Note that using custom cutoff points is
%       about n times slower than equal spacing. That means that, compared
%       to the equally-spaced levels case, where the time doesn't scale
%       with n, when custom cutoffs are used the time is O(n).

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('cutoffs', [], @(x) isempty(x) || (isvector(x) && isnumeric(x)));

% show defaults if requested
if nargin == 1 && strcmp(img, 'defaults')
    parser.parse;
    disp(parser.Results);
    clear img;
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

if isempty(params.cutoffs)
    img = min(max(floor(n*img)/(n-1), 0), 1);
else
    imgOut = zeros(size(img));
    % everything smaller than the first cutoff should be set to 0
    for i = 2:length(params.cutoffs)
        imgOut(img >= params.cutoffs(i-1) & img < params.cutoffs(i)) = (i-1)/(n-1);
    end
    % everything larger than the last cutoff should be set to 1
    imgOut(img >= params.cutoffs(end)) = 1;
    
    img = imgOut;
end

if ismatrix(img)
    crop = [1 1 size(img)];
else
    crop = [];
end

end
