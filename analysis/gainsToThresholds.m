function thresholds = gainsToThresholds(L, directions)
% gainsToThresholds Calculate thresholds in some directions given the gain
% matrix.
%   thresholds = gainsToThresholds(L, directions) calculates the threshold
%   magnitudes in the given `directions` (a cell array of vectors or a
%   single such vector), assuming that, after being transformed by the gain
%   matrix `L`, the threshold surface is a hypersphere of unit radius
%   centered at the origin.

if ~iscell(directions)
    directions = {directions};
end

% threshold location: t*v, where v is direction vector (not necessarily
% normalized)
% transformed threshold location: t*L*v
% hypersphere: t^2*v'*L'*L*v = 1
% so t = 1 / sqrt(v'*(L'*L)*v) = 1/norm(L*v);

thresholds = 1./ cellfun(@(v) norm(L*v(:)), directions);

end