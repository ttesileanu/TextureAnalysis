function eivec=flfs_findeivec(configs,subcat_freqs,flfs)
% function eivec=flfs_findeivec(configs,subcat_freqs,flfs): computes stable eigenvector of a propagation matrix
%
% configs: output of flfs_enumerate, listing all configurations of a
%    falling-leaves/sticks setup
% subcat_freqs:  a vector of relative frequencies of the subcategories, need not add to 1
% flfs: structure defining falling-leaves/sticks categories and subcategories, returned by flfs_define
%
% eivec: column vector of length length(configs); length sums to 1
%
% strategy:  first find eigenvector of the matrix of permutation
% transitions, and then sum over permutations that are all parents of the
% same configuration
%
% See also:  FLFS_SEQREDUCE, FLFS_DEFINE, FLFS_ENUMERATE, FLFS_MAKEPROPAG, FLFS_ENUMERATE_DEMO, FLFS_EIVEC_TEST.
%
nconfigs=length(configs);
%
eivec=zeros(nconfigs,1);
for iconfig=1:nconfigs
    parents=configs(iconfig).parents;
    for ipar=1:size(parents,1)
        freqs=subcat_freqs(parents(ipar,:)); % the frequencies in the permutation, most recent first
        unused_freqs=sum(freqs)-cumsum([0 freqs(1:end-1)]);
        nonz=find(unused_freqs>0);
        prob_perm=prod(freqs(nonz))/prod(unused_freqs(nonz)); %this is the probability of the permutation
        eivec(iconfig)=eivec(iconfig)+prob_perm;
    end
end
eivec=eivec/sum(eivec);
return
