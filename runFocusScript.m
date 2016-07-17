function dataNI = runFocusScript(dataNI, varargin)
% runFocusScript Fit mixture of Gaussians to image patch statistics, find
% in-focus patches.
%   dataNI = runFocusScript(dataNI) fits a mixture of Gaussian to the image
%   patches in each of the structures from the array dataNI.indA. If the
%   'focusImage' options is used (see below), the function also identifies
%   in-focus patches.
%
%   The Gaussian mixture object is stored in a field 'obj' of the output
%   dataNI.indA(i) structure, while the cluster to which each patch is most
%   likely to belong is stored in the vector 'cx' and the Mahalanobis
%   distance to each Gaussian is stored in the matrix 'mahal'. The mean of
%   the cluster (cx) values for all of the image patches for every image is
%   stored in 'mn'.
%
%   When the 'focusImage' option is used (see below), a member 'focus' is
%   added to each of the entries in the output dataNI.indA structure array.
%   It is a structure with members
%    chist:
%       A vector describing, for each patch of the focus image, which of
%       the Gaussian mixture components it belongs to.
%    component:
%       The Gaussian mixture component that is judged to be in-focus. This
%       is the one that has the majority of entries in the chist vector.
%    ev:
%       Subset of the 'ev' matrix corresponding to the in-focus patches.
%       (see analyzeImageSetModNoPC)
%    covM:
%       Covariance matrix for the ev matrix restricted to the in-focus
%       patches.
%
%   Options:
%    'focusImage': int
%       The index of an image that has most of its patches in focus. This
%       is used to identify which of the two Gaussians clusters is most
%       likely to contain the in-focus patches.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('focusImage', [], @(n) isscalar(n) && isnumeric(n));

% parse
parser.parse(varargin{:});
params = parser.Results;

% run 2 Gaussian decomposition on all analyses
for i=1:length(dataNI.indA)
    disp(['Running Focus Analysis ',num2str(i)])
    options = statset('Display','off');
    dataNI.indA(i).obj = gmdistribution.fit(dataNI.indA(i).ev,2, ...
        'options',options,'SharedCov',false,'replicates',10);
    
    % assign every image patch to one of the two gaussians
    dataNI.indA(i).cx = dataNI.indA(i).obj.cluster(dataNI.indA(i).ev);
    dataNI.indA(i).mahal = dataNI.indA(i).obj.mahal(dataNI.indA(i).ev);
    
    % for every image, compute how many patches are in cluster 1 or cluster 2
    % on average
    for j=1:numel(dataNI.indA(i).ic.name);
        ix=(dataNI.indA(i).ic.image==j);
        dataNI.indA(i).mn(j)=mean(dataNI.indA(i).cx(ix));
    end
    
    if ~isempty(params.focusImage)
        cs = dataNI.indA(i).cx(dataNI.indA(i).ic.image == params.focusImage);
        dataNI.indA(i).focus.chist = cs;
        if isempty(cs)
            warning([mfilename ':nofocus'], ['Focus image index ' int2str(params.focusImage) ...
                ' does not have any patches in analysis ' int2str(i) '.']);
        else
            c = round(dataNI.indA(i).mn(params.focusImage));
            
            disp(['Using image ' int2str(params.focusImage) ', ' dataNI.indA(i).ic.name{params.focusImage} ...
                ' as in-focus for analysis ' int2str(i) '.']);
            
            dataNI.indA(i).focus.component = c;
            dataNI.indA(i).focus.ev = dataNI.indA(i).ev(dataNI.indA(i).cx==c,:);
            dataNI.indA(i).focus.covM = squeeze(dataNI.indA(i).obj.Sigma(:,:,c));
        end
    end
end

end