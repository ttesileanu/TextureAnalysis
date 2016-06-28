% generateStatFiles Generate texture statistics file.
%   This will load the statistics files if they already exist, unless
%   params.forcestats is true. If params.forcestats is true or the
%   statistics files cannot be found or loaded, they will be generated.

Ffilter = cell(numel(patchSizes)*numel(blockAFs));

crt = 1;
for i=1:numel(patchSizes)
    for j=1:numel(blockAFs)
        outputTag = strcat('_N0',num2str(blockAFs(j)),...
            '_PCfrac',num2str(100*PCfrac),'_PS',num2str(patchSizes(i)),'.mat');        
        outputFname = fullfile(analysesDirectory, analysesFileName, outputTag);
        if ~params.forcestats
            try
                rawdata = open(outputFname);
                data = rawdata.data;
                Ffilter{crt} = rawdata.crtFfilter;
                disp(['Loaded from file Patch Size ',num2str(patchSizes(i)),' , Block Avg Factor ', num2str(blockAFs(j))]);
            catch
                data = [];
            end
        else
            data = [];
        end
        
        if isempty(data)
            disp(['Running Patch Size ',num2str(patchSizes(i)),' , Block Avg Factor ', num2str(blockAFs(j))]);
            data=struct();
            [ covM, ev, imgXY, Ffilter{crt}, sharpness ] = ...
                analyzeImageSetModNoPC(imgNamesFile, imgDirectory,...
                    imgNamesFile_filter, imgDirectory_filter, blockAFs(j), ...
                    PCfrac, patchSizes(i));
            data.covM = covM;
            data.ev = ev;
            data.imageCoordinates = imgXY;
            data.PCfrac = PCfrac;
            data.patchSize = patchSizes(i);
            data.blockAvgFactor = blockAFs(j);
            data.sharpness = sharpness;
            
            crtFfilter = Ffilter{crt}; %#ok<NASGU>
            save(outputFname, {'data', 'crtFfilter'});
        end
        crt = crt + 1;
    end
end

clear blockAFs covM data ev i imgXY j patchSizes PCfrac outputTag crtFfilter crt rawdata