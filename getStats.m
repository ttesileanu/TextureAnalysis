function [ P, ev, entropy ] = getStats( imgs )
% getStats Calculate statistics from image patches.
%   [P, ev, entropy] = getStats(imgs) calculates statistics for the image
%   patches given in the 3-dimensional array imgs. The last index indexes
%   the patches.
%
%   Returns:
%    P:
%       Probabilities for each of the 16 glider configurations calculated
%       for each image patch. (see processBlock)
%    ev:
%       Values of the 10-dimensional independent combinations of
%       probabilities for each image patch. (see getStatistics and
%       processBlock).
%    entropy:
%       Entropy values for each of the image patches. (see processBlock)
%
%   See also: processBlock, getStatistics.

P=zeros(size(imgs,3),16);
ev=zeros(size(imgs,3),10);
entropy=zeros(size(imgs,3),1);
for i=1:size(imgs,3)
    [P(i,:), ev(i,:), entropy(i)]=processBlock(imgs(:,:,i));
end

end

