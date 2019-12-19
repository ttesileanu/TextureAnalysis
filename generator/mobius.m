function m=mobius(n)
% m=mobius(n) is the mobius function
% mobius(1)=1
% if n is a product of k distinct primes, mobius(n)=(-1)^k
% otherwise m=0
%
%   See also:  TOTIENT.
%
m=0;
if n<1; return; end
if n>floor(n); return; end
if n==1; m=1; return; end
f=factor(n);
uf=unique(f);
if n>prod(uf); return; end %n has  a square factor
m=(-1)^length(uf);
return
end
