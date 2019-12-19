function h=histinfo(pvec)
%
% function h=histinfo(pvec) gives the histogram information, in bits
% sum(pvec) assumed to be 1. but needn't be 1-dimensional. zero and negative elements ignored, 
% lacks all the bells and whistles of histent.
%
% if all elements are known positive, then use histinfo_nz, which is faster.
% See also TBLXINFO, HISTENT, HISTINFO_NZ.
%
pnz=pvec(:);
pnz=pnz(find(pnz>0)); %revised 7 Jul 17 to omit the test for empty pnz; an empty pnz will still return 0
%if (isempty(pnz))
%   h=0;
%else
%   h=-pnz*log(pnz)'/log(2);
%end
h=-pnz'*log(pnz)/log(2);
return

