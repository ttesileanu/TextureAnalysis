function dataNI = convertData(old_dataNI)
% convertData Change some field names in dataNI and exchange some of the
% entries of the covariance matrix to match the conventions from the
% psychophysics data.
%   dataNI = convertData(old_dataNI) slightly changes the old_dataNI
%   structure to use naming conventions and column/row orderings that match
%   previous work.
%
%   'dataNI' will have a single member, a structure array called
%   'indA'. The fields of this structure are:
%    ev:
%       Values of the 10 independent texture parameters for all the image
%       patches. (see analyzeImageSetModNoPC)
%    covM:
%       Covariance matrix for the texture parameters.
%       NOTE: compared to analyzeImagesSetModNoPC, indices 7 and 8 are
%             flipped in this script's covM matrices. This matches the
%             order used in the psychophysics experiments.
%    ic:
%       Contains information regarding the images from which the patches
%       come, and the location of each patch within the image.
%    N:
%       Block average factor used in the analysis.
%    R:
%       Patch size (after block averaging).
%    sharpness:
%       Vector of sharpness values for each image patch. (see getImgStats)
%
%   See also: generateStatFiles, analyzeImageSetModNoPC, getImgStats.

dataNI = struct;
for i = 1:length(old_dataNI.indA)
    old_data = old_dataNI.indA(i);
    
    new_data = rmfield(old_data, {'imageCoordinates', 'blockAvgFactor', 'patchSize'});
    
    new_data.ic     = old_data.imageCoordinates;
    new_data.N      = old_data.blockAvgFactor;
    new_data.R      = old_data.patchSize;
    
    %reorder covM to align with psychophysics
    new_data.covM(:,[7 8])  = new_data.covM(:,[8 7]);
    new_data.covM([7 8],:)  = new_data.covM([8 7],:);
    
    dataNI.indA(i) = new_data;
end

end