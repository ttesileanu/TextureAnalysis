function u=outprod(v,w)
%
% u=outprod(v,w) computes the outer product of v and w
%
% size(u)=[sizetrim(v) sizetrim(w)] where
%    sizetrim(v)=size(v) with trailing 1's deleted
%
% typically v,w are column vectors, but either can be multidimensional
%
%   See also:  SIZETRIM, OUTPROD_CONTRACT, KRON.
%
vq=repmat(v,[ones(1,length(sizetrim(v))),sizetrim(w)]);
wr=reshape(w,[ones(1,length(sizetrim(v))),sizetrim(w)]);
wq=repmat(wr,[sizetrim(v),ones(1,length(sizetrim(v)))]);
u=vq.*wq;
return
