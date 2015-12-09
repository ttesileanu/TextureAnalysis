%generateStatFiles

PCfrac = 1;
patchSizes = [32,48];
%patchSizes = [32];
%blockAFs = [2,4];
blockAFs = [2];

Ffilter = cell(numel(patchSizes)*numel(blockAFs));

crt = 1;
for i=1:numel(patchSizes)
    for j=1:numel(blockAFs)
        disp(['Running Patch Size ',num2str(patchSizes(i)),' , Block Avg Factor ', num2str(blockAFs(j))])
        data=struct();
        [ covM, ev, imgXY, Ffilter{crt}, sharpness ] = analyzeImageSetModNoPC( imgNamesFile, imgDirectory,...
            imgNamesFile_filter, imgDirectory_filter, blockAFs(j), PCfrac, patchSizes(i));
        data.covM = covM;
        data.ev  =ev;
        data.imageCoordinates = imgXY;
        data.PCfrac = PCfrac;
        data.patchSize = patchSizes(i);
        data.blockAvgFactor = blockAFs(j);
        data.sharpness = sharpness;
        
        outputTag = strcat('_N0',num2str(blockAFs(j)),...
            '_PCfrac',num2str(100*PCfrac),'_PS',num2str(patchSizes(i)),'.mat');
        save(strcat(analysesDirectory,analysesFileName,outputTag),'data');
        
        crt = crt + 1;
        
    end
    
end

clear blockAFs covM data ev i imgXY j patchSizes PCfrac outputTag