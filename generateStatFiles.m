function dataNI = generateStatFiles(params)
% generateStatFiles Generate texture statistics file.
%   dataNI = generateStatFiles(params) calculates the statistics given the
%   parameters that were passed to a processData script. This function also
%   stores the statistics to file, provided params.subSelect is empty.
%
%   See also: processData.

if isempty(params.filters) || strcmp(params.filters, 'load')
    Ffilter = cell(numel(params.patchSizes)*numel(params.blockAFs));
else
    Ffilter = params.filters;
end

dataNI = struct;

crt = 1;
for i=1:numel(params.patchSizes)
    patchSize = params.patchSizes(i);
    for j=1:numel(params.blockAFs)
        blockAF = params.blockAFs(j);
        outputTag = strcat('_N0', int2str(blockAF), ...
            '_PS', int2str(patchSize), '.mat');
        outputFname = fullfile(params.analysesDirectory, ...
            strcat(params.outPrefix, outputTag));
        
        % does the file already exist?
        try
            rawData = open(outputFname);
        catch
            rawData = [];
        end

        if ~isempty(rawData)
            if strcmp(params.filters, 'load')
                Ffilter{crt} = rawData.data.Ffilter;
                disp(['Loading filter from ' outputFname '.']);
            end
        end
        
        doGenerate = true;
        if ~isempty(rawData) && isempty(params.subSelect) && params.loadStats
            canLoad = true;
            % perhaps we can load from file?
            try
                if ~isempty(Ffilter{crt}) && ~isequal(Ffilter{crt}, rawData.data.Ffilter)
                    canLoad = false;
                else
                    if ~isequal(rawData.data.patchSize, patchSize) || ...
                       ~isequal(rawData.data.blockAvgFactor, blockAF)
                        canLoad = false;
                    end
                end
            catch
                canLoad = false;
            end
            
            if canLoad
                data = rawData.data;
                doGenerate = false;
                disp(['Loaded statistics for Patch Size ', int2str(patchSize), ', Block AF ', int2str(blockAF) ...
                    ' from ' outputFname '.']);
            end
        end
        
        if doGenerate
            disp(['Running Patch Size ', int2str(patchSize), ', Block AF ', int2str(blockAF)]);
            
            segmentation = struct(...
                'segs', params.segs, ...
                'segSel', params.segSel, ...
                'segField', params.segField, ...
                'segFct', params.segFct);
            
            analysisOpts = {...
                'segmentation', segmentation, ...
                'images', params.images, ...
                'subSelect', params.subSelect, ...
                'filterFull', params.filterFull, ...
                'overlappingPatches', params.overlappingPatches, ...
                'fullImageEv', params.fullImageEv, ...
                'progressStart', params.progressStart, ...
                'progressEvery', params.progressEvery};
            if ~isempty(Ffilter{crt})
                analysisOpts = [analysisOpts {'filter' Ffilter{crt}}]; %#ok<AGROW>
            end
            [ covM, ev, others ] = ...
                analyzeImageSetModNoPC(params.imgNamesFile, params.imgDirectory,...
                params.imgNamesFile_filter, params.imgDirectory_filter, blockAF, ...
                patchSize, analysisOpts{:});
            
            data = struct;
            
            data.covM = covM;
            data.ev = ev;
            data.patchSize = patchSize;
            data.blockAvgFactor = blockAF;
            
            data.imageCoordinates = others.imageCoordinates;
            data.sharpness = others.sharpness;
            data.patchStats = others.patchStats;
            
            extraFields = {'evF', 'evB', 'covF', 'covB', ...
                'origImages', 'blockAFImages', 'whitenedImages', 'binarizedImages', ...
                'segmentationMasks'};
            for k = 1:length(extraFields)
                if isfield(others, extraFields{k})
                    data.(extraFields{k}) = others.(extraFields{k});
                end
            end
            
            data.params = params;
            if isfield(data.params, 'segs')
                % don't store the whole segmentation data in the output
                data.params = rmfield(data.params, 'segs');
            end
            if isfield(data.params, 'filters')
                % storing all the filters would take too much space
                data.params = rmfield(data.params, 'filters');
            end
            
            Ffilter{crt} = others.Ffilter;
            data.Ffilter = Ffilter{crt};
            
            if (isempty(rawData) || params.forceSave) && isempty(params.subSelect)
                save(outputFname, 'data');
                disp(['Saved statistics for Patch Size ', int2str(patchSize), ', Block AF ', int2str(blockAF) ...
                    ' to ' outputFname '.']);
            else
                if ~isempty(rawData) && isempty(params.subSelect)
                    warning([mfilename ':nooverw'], ['Can''t save analysis for patch size ' int2str(patchSize) ...
                        ', block AF ' int2str(blockAF)  ' -- file already exists.']);
                end
            end
        end
        
        dataNI.indA(crt) = data;
        crt = crt + 1;
    end
end

end