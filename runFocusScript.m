% runFocusScript Fit mixture of Gaussians to image patch statistics.
%   The script fits a mixture of two Gaussians to the image patches in each
%   of the structures from the array dataNI.indA. The gaussian mixture
%   object is stored in a field 'obj' of dataNI.indA(i), while the cluster
%   to which each patch is most likely to belong is stored in the vector
%   'cx' and the Mahalanobis distance to each Gaussian is stored in the
%   matrix 'mahal'. The mean of the cluster (cx) values for all of the
%   image patches for every image is stored in 'mn'.

% run 2 Gaussian decomposition on all analyses

Nimg = numel(unique(dataNI.indA(1).ic.image));

for i=1:length(dataNI.indA)
    disp(['Running Focus Analysis ',num2str(i)])
    options = statset('Display','off');
    dataNI.indA(i).obj = gmdistribution.fit(dataNI.indA(i).ev,2,'options',options,'SharedCov',false,'replicates',10);
    
    % assign every image patch to one of the two gaussians
    dataNI.indA(i).cx = dataNI.indA(i).obj.cluster(dataNI.indA(i).ev);
    dataNI.indA(i).mahal = dataNI.indA(i).obj.mahal(dataNI.indA(i).ev);
    
    % for every image, compute how many patches are in cluster 1 or cluster 2
    % on average
    for j=1:Nimg
        ix=(dataNI.indA(i).ic.image==j); 
        dataNI.indA(i).mn(j)=mean(dataNI.indA(i).cx(ix));
    end
end

clear i ix j Nimg options