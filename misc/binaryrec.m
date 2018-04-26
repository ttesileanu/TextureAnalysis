function vecs = binaryrec(mags, uvecs)
% BINARYREC Recompose magnitude+direction pairs to vectors in binary
% texture groups.
%   vecs = BINARYREC(mags, uvecs) recomposes vectors decomposed into
%   magnitude+direction by BINARYDEC. The vectors are scaled such that the
%   direction from the [1/2, 1/2, ..., 1/2] origin of the coordinate system
%   to the vector is preserved.
%
%   The directions `uvecs` can be given either as an n x 2k matrix, or as a
%   cell array of 2k-dimensional vectors. In either case, the output is
%   returned as an n x 2k matrix.
%
%   See also: BINARYDEC.

% handle cell-array input
if iscell(uvecs)
    uvecs = cell2mat(cellfun(@(v) v(:)', uvecs(:), 'uniform', false));
end

vecs = 1/2 + bsxfun(@times, uvecs - 1/2, mags(:));

end