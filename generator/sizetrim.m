function s=sizetrim(v)
%
% s=sizetrim(v) returns size(v) with trailing 1's deleted
%   except that sizetrim(singleton)=1
%
%   See also:  OUTPROD, SQUEEZE.
%
s=size(v);
k=find(s~=1);
if (length(k)==0)
   s=1;
   return
end
s=s(1:max(k));
return
