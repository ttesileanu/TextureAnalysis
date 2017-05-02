function [segs, segSel] = groupPennSegs(image_names, seg_path)

segSel = false(length(image_names), 1);
segs = struct;
for i = 1:length(image_names)
    crt_name = image_names{i};
    [crt_path, crt_file] = fileparts(crt_name);
    % remove '_LUM' from image file names, as it doesn't appear in segs
    if length(crt_file) > 4 && strcmp(crt_file(end-3:end), '_LUM')
        crt_seg_file0 = crt_file(1:end-4);
    end
    % sometimes extra label were added to the name of the segmented file
    % use 'dir' to find the right segmentation
    crt_seg_file = [crt_seg_file0 '_segmented'];
    possible = dir(fullfile(seg_path, crt_path, [crt_seg_file '*.mat']));
    % exclude directories
    possible = possible(arrayfun(@(s) ~s.isdir, possible));
    segs(i).image = crt_name;
    if ~isempty(possible)
        % this image does have a segmentation
        if length(possible) > 1
            % maybe too many...
            warning(['Image ' crt_name ' has several segmentations. Using ' ...
                fullfile(seg_path, crt_path, possible(1).name) '.']);
        end
        crt_seg_fullfile = fullfile(seg_path, crt_path, possible(1).name);
        crt_seg_data = open(crt_seg_fullfile);
        
        crtFG.segPath = fullfile(crt_path, possible(1).name);
        crtFG.objMat = crt_seg_data.segmentation;
        crtFG.objMatLabels = crt_seg_data.tags;
        
        segSel(i) = true;
    else
        crtFG.segPath = '';
        crtFG.objMat = [];
        crtFG.objMatLabels = {};
    end
    segs(i).FG = crtFG;
end

end