% load image analyses
dr0 = analysesDirectory;
dr  = dir(strcat(analysesDirectory,analysesFileName,'*.mat'));

dataNI=struct;
disp(['Number of analyses: ' num2str(length(dr))])
for i=1:length(dr),
    load(strcat(dr0,dr(i).name));
    dataNI.indA(i).ev   = data.ev;
    
    %reorder covM to align with psychophysics
    data.covM(:,[7 8])  = data.covM(:,[8 7]);
    data.covM([7 8],:)  = data.covM([8 7],:);
    dataNI.indA(i).covM = data.covM;
    
    dataNI.indA(i).ic   = data.imageCoordinates;
    dataNI.indA(i).N    = data.blockAvgFactor;
    dataNI.indA(i).R    = data.patchSize;
    dataNI.indA(i).sharpness = data.sharpness;
    disp(['Loading Analysis ',num2str(i)])
end

clear data dr dr0 i