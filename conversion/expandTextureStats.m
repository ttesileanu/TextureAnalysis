function ev = expandTextureStats(ev, G, type)
% expandTextureStats Expand grayscale texture statistics from an
% independent form to a probability form.
%   evP = expandTextureStats(ev, G) expands the N x S matrix of `G`-level
%   grayscale statistics `ev` by adding the missing probability values such
%   that each consecutive set of `G` values adds up to 1.
%
%   For `G == 2`, this makes sure to disentangle the differences calculated
%   in that case (see `analyzeTexture`).
%
%   expandTextureStats(ev, G, '3d') returns the stats as an N x G x nPlanes
%   3d array (i.e., without flattening the last two dimensions).

if nargin < 3
    type = '2d';
end

% handle the `G == 2` case
if G == 2
    % alpha and theta have minus signs for historical reasons
    ev(:, [1 6:9]) = -ev(:, [1 6:9]);
    % need to replace differences (nEven - nOdd)/nTotal by just first
    % components, nEven/nTotal.
    ev = (1 + ev) / 2;
    % also need to exchange some entries!
    reorder = [1 3 2 5 4 6 8 9 7 10];
    ev = ev(:, reorder);
end

% easiest way is to go through a 3-dimensional array
nPlanes = size(ev, 2)/(G-1);
ev = reshape(ev, [], G-1, nPlanes);
ev(:, G, :) = 1 - sum(ev, 2);

switch type
    case '2d'
        ev = reshape(ev, [], nPlanes*G);
    case '3d'
    otherwise
        error([mfilename ':badtype', 'Unrecognized type argument.');
end

end