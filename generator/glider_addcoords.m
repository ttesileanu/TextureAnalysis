function [g,optsused]=glider_addcoords(gin,opts)
% g=glider_addcoords(gin) adds the indexing coordinates to a glider structure
%
% gin:  a glider structure
% gin.matrix: array of 0's and 1's
% opts: options
%   (none used at present)
%
% g: has all the fields of gin and also
% g.size:  size of bounding box (rows and cols)
% g.rank:  number of occupied cells
% g.inds:  pointers (1-based) into occupied cells; size(g.inds)=[g.rank 2]
%          sorted by rows
%
% optsused: options used
%
%   1-d gliders were transposed prior to 12/28/06, fixed
%   See also:  GLIDER_TEST, GLIDER_CHOOSE.
%
if (nargin<=1) opts=[]; end
optsused=opts;
%
g=gin;
g.size=size(gin.matrix);
g.rank=sum(gin.matrix(:));
if (size(gin.matrix,2)==1)  % was (gin.matrix,1) until 12/28/06
    rows=reshape(find(gin.matrix),length(find(gin.matrix)),1);
    cols=ones(length(rows),1);
elseif (size(gin.matrix,1)==1) % was (gin.matrix,2) until 12/28/06
    cols=reshape(find(gin.matrix),length(find(gin.matrix)),1);
    rows=ones(length(cols),1);
else
    [rows,cols]=find(gin.matrix);
end
g.inds=sortrows([rows cols]);
%
return
