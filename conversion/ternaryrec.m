function vecs = ternaryrec(mags, uvecs)
% TERNARYREC Recompose magnitude+direction pairs to vectors in ternary
% texture planes.
%   vecs = TERNARYREC(mags, uvecs) recomposes vectors decomposed into
%   magnitude+direction by TERNARYDEC. The vectors are scaled such that the
%   direction from the [1/3, 1/3, ..., 1/3] origin of the coordinate system
%   to the vector is preserved.
%
%   The directions `uvecs` can be given either as an n x 3k matrix, or as a
%   cell array of 3k-dimensional vectors. In either case, the output is
%   returned as an n x 3k matrix.
%
%   See also: TERNARYDEC.

% handle cell-array input
if iscell(uvecs)
    uvecs = cell2mat(cellfun(@(v) v(:)', uvecs(:), 'uniform', false));
end

vecs = 1/3 + bsxfun(@times, uvecs - 1/3, mags(:));

end