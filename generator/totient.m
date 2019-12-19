function t=totient(d,b)
% t=totient(d,b) computes the generalized totient function
% for b=0 or omitted, t=totient(d)
% note also that totient(d,1)=mobius(d)
%
% t=d*sum(f|d, bf==0 mod d)*mobius(f)*(1/f)
%
% See ...\stim\graylevs_mtc.doc for T'(D,b), a use of the generalized totient function
%
%   See also:  MOBIUS.
%
if (nargin<=1)
    b=0;
end
t=0;
if d<1; return; end
if d>floor(d); return; end
if d==1; t=1; return; end
f1=find(mod(d,[1:d])==0); %the divisors of d
f2=find(mod(b*[1:d],d)==0); %bf==0 (mod d)
flist=intersect(f1,f2);
t=0;
for pf=1:length(flist)
    f=flist(pf);
    t=t+(d/f)*mobius(f);
end
return
end

