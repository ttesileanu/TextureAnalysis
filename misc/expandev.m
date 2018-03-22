function ev = expandev(ev, G)
% EXPANDEV Expand grayscale texture statistics from an independent form to
% a probability form.
%   evP = EXPANDEV(ev, G) expands the N x S matrix of `G`-level grayscale
%   statistics `ev` by adding the missing probability values such that each
%   consecutive set of `G` values adds up to 1.

% easiest way is to go through a 3-dimensional array
n_planes = size(ev, 2)/(G-1);
ev = reshape(ev, [], G-1, n_planes);
ev(:, G, :) = 1 - sum(ev, 2);
ev = reshape(ev, [], n_planes*G);

end