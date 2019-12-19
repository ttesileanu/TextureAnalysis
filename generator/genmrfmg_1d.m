function rowvec=genmrfmg_1d(n,pixelprobs,dyadprobs,firstpxl)
% rowvec=genmrfmg_1d(n,pixelprobs,dyadprobs,firstpxl)
% make a 1d Markov process with pixelprobs for single-pixel probabilities
% and dyadprobs for the dyad probabilities
%   pixelprobs: size is [1 ng]
%   dyadprobs: size is [ng ng]
%   n: number of checks to generate
%   firstpxl:  if supplied and not empty, the value of the first pixel (added 21-Sep-12)
%      otherwise this is generated according to the stable probabilities
%
%    See also:  GENMRFMG, GENMRFM_1D.
%
ng=length(pixelprobs);
if nargin<=3
    firstpxl=[];
end
for k=1:n
    if (k==1)
        if (isempty(firstpxl))
            rowvec(1,k)=sum(rand(1)>cumsum(pixelprobs));
        else
            rowvec(1,k)=firstpxl;
        end
    else
        probs=dyadprobs(rowvec(k-1)+1,:);
        if (sum(probs)==0)
            rowvec(1,k)=sum(rand(1)>cumsum(pixelprobs));
        else
            rowvec(1,k)=sum(rand(1)>(cumsum(probs)/sum(probs)));
        end
    end
    rowvec(1,k)=min(rowvec(1,k),ng-1);
end
return
