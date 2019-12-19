function unbiased_list=mtc_unbiased(ng)
% unbiased_list=mtc_unbiased(ng) is a utility that generates the unbiased
% vectors for each base harmonic
%
% unbiased_list(baseharm,:): the unbiased vectors for baseahrm
%   dim 1 of unbiased_list runs up to the largest factor of ng
%   dim 2 of unbiased_list is length ng
% for ng prime, it is ones(1,ng)/ng
%
if (isprime(ng))
    unbiased_list=ones(1,ng)/ng;
else
    for k=1:ng/2
        if mod(ng,k)==0
            unbiased_list(k,:)=(double(mod([0:ng-1],k)==0))/(ng/k);
        end
    end
end
return

