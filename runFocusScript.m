% run 2 Gaussian decomposition on all analyses

Nimg = numel(unique(dataNI.indA(1).ic.image));

for i=1:length(dataNI.indA)
    disp(['Running Focus Analysis ',num2str(i)])
    options = statset('Display','off');
    dataNI.indA(i).obj = gmdistribution.fit(dataNI.indA(i).ev,2,'options',options,'SharedCov',false,'replicates',10);
    
    % assign every image patch to one of the two gaussians
    dataNI.indA(i).mahal = zeros(size(dataNI.indA(i).ev, 1), 2);
    for j=1:size(dataNI.indA(i).ev,1)
        dataNI.indA(i).cx(j)=dataNI.indA(i).obj.cluster(dataNI.indA(i).ev(j,:));
        dataNI.indA(i).mahal(j, :) = dataNI.indA(i).obj.mahal(dataNI.indA(i).ev(j, :));
    end
    
    % for every image, compute how many patches are in cluster 1 or cluster 2
    % on average
    for j=1:Nimg
        ix=(dataNI.indA(i).ic.image==j); 
        dataNI.indA(i).mn(j)=mean(dataNI.indA(i).cx(ix));
    end
end

clear i ix j Nimg options