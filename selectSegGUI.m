function segSel = selectSegGUI(segData, imgDirectory, varargin)
% selectSegGUI Select best segmentation for each image, and choose images
% to discard because of bad segmentations.
%   segSel = selectSegGUI({segs, segSel, imgNames}, imgDirectory) opens the
%   images whose names are given in the imgNames cell array, searching for
%   them in the imgDirectory folder, and overlays the segmentations from
%   the segs structure. Each image can have several possible segmentations,
%   and the segSel vector indicates which segmentation to use for each
%   image; a NaN value indicates an image that should be skipped. The GUI
%   allows changing the selections in the segSel vector. After the GUI ends,
%   the final segSel vector is returned by the function. Set the segSel
%   input argument to an empty vector to default to using the first
%   segmentation for each image.
%
%   The segs variable contains the segmentation data as a cell array. Each
%   element of this cell array is structure array, each entry of which is a
%   possible segmentation. The segs variable can also contain several
%   different types of segmentation, such as foreground, per-object
%   segmentations, or contour data. These are identified by different
%   fields in the structure array, and can be chosen using a dropdown box
%   in the GUI. All the fields of the segs{1} structure array that are
%   numeric matrices are considered to be valid segmentation choices.
%
%   segSel = selectSegGUI({segs, segSel}, imgDirectory) uses a different
%   representation of the segmentations. In this case, segs is a structure
%   array with as many entries as there are images, and it has two fields.
%   One is called 'image', 'images', 'name', or 'names', and it provides
%   the names of the images that the segmentation refers to (like imgNames
%   above). The other field's name is unimportant; it provides the same
%   data as each entry of the segs cell array described above, but in the
%   form of a structure array.
%
%   selectSegGUI(segs, imgDirectory) works as above, but starts with the
%   default segSel (all ones).
%
%   segSel = selectSegGUI(segFileName, imgDirectory) loads the segmentation
%   data from the given file. The file should contain at least the
%   segmentation data in a variable called 'segs' or 'segmentations'. This
%   can be in the structure array format described above. If the
%   segmentation data is in cell array format instead, then the image names
%   must be provided separately, in a variable called 'names',
%   'image_names', or 'imageNames'. Finally, the initial segmentation
%   choices are loaded from a variable called 'segSel', 'seg_select',
%   'segSelect', 'seg_selections', or 'segSelections', if available.

% load previous settings
try
    settings = open('segGuiSettings.mat');
    settings = settings.settings;
catch
    settings = struct;
end

settings = setIfMissing(settings, ...
    'figPos', [0.1 0.05 0.85 0.85], ...
    'singleImagePos', [0.2, 0.1, 0.6, 0.7], ...
    'field', 'fgMat', ...
    'showSeg', true, ...
    'saveChanges', true, ...
    'confirmSave', true, ...
    'minRowSpace', 1, ...   % in character heights (see Matlab documentation for Units)
    'minColSpace', 2, ...   % in character widths
    'rowHeight', 6, ...     % in character heights
    'maxAspect', 2, ...     % horizontal size of thumbnails is capped at maxAspect * row height
    'titleHeight', 1, ...   % height reserved for image title (in character heights)
    'controlsHeight', 1, ...% space reserved for picture controls (in character heights)
    'figEdge', 1, ...       % figure edge in character widths
    'topEdge', 2, ...       % extra edge at top of figure (in character heights)
    'bottomEdge', 2, ...    % extra edge at bottom of figure (in character heights)
    'leftEdge', 0, ...      % extra edge at figure left (in character widths)
    'rightEdge', 0 ...      % extra edge at figure right (in character widths)
);

% override settings given the extra arguments
if mod(numel(varargin), 2) ~= 0
    error([mfilename ':badargs'], 'Additional arguments should always come in key/value pairs.');
end
for i = 1:numel(varargin)/2
    varKey = varargin{2*i-1};
    varValue = varargin{2*i};
    
    if ~ischar(varKey) || ~isvector(varKey)
        error([mfilename ':badkey'], 'Argument names should be strings.');
    end
    
    if ~isfield(settings, varKey)
        error([mfilename ':badarg'], ['Unknown key ''' varKey '''.']);
    end
    settings.(varKey) = varValue;
end

% handle the various input formats
% XXX maybe I should move some of this logic to a function
fromFile = false;
if ischar(segData)
    % load from file
    fromFile = true;
    segFileName = segData;
    segFileData = open(segData);
    
    possibleSegFields = {'segs', 'segmentations'};
    fieldIdx = find(isfield(segFileData, possibleSegFields), 1);
    if isempty(fieldIdx)
        error([mfilename ':fnoseg'], 'No segmentation data found in file.');
    end
    segFieldNameInFile = possibleSegFields{fieldIdx};
    prelimSegs = segFileData.(segFieldNameInFile);
    
    if iscell(prelimSegs)
        % need separate image name data
        possibleImageFields = {'names', 'image_names', 'imageNames'};
        fieldIdx = find(isfield(segFileData, possibleImageFields), 1);
        if isempty(fieldIdx)
            error([mfilename ':fnoimgs'], 'Can''t find image names in file.');
        end
        imageField = possibleImageFields{fieldIdx};
        prelimImgNames = {segFileData.(imageField)};
    else
        prelimImgNames = [];
    end
    
    % check if we have any segSel data
    possibleSegSelFields = {'segSel', 'seg_select', 'segSelect', ...
        'seg_selections', 'segSelections'};
    fieldIdx = find(isfield(segFileData, possibleSegSelFields), 1);
    if isempty(fieldIdx)
        prelimSegSel = [];
        segSelFieldNameFile = possibleSegSelFields{1};
    else
        segSelFieldNameFile = possibleSegSelFields{fieldIdx};
        prelimSegSel = segFileData.(segSelFieldNameFile);
    end
    
    % format this so that it's processed by the code below
    if ~isempty(prelimImgNames)
        segData = {prelimSegs, prelimSegSel, prelimImgNames};
    else
        segData = {prelimSegs, prelimSegSel};
    end
elseif ~iscell(segData)
    segData = {segData};
end
if iscell(segData)
    if numel(segData) == 3
        % this is the cell array format
        segs = segData{1};
        segSel = segData{2};
        imgNames = segData{3};
    else
        % segmentation and image names combined in one structure
        possibleImageFields = {'image', 'images', 'name', 'names'};
        fieldIdx = find(isfield(segData{1}, possibleImageFields), 1);
        if isempty(fieldIdx)
            error([mfilename ':nonames'], 'Can''t find image names.');
        end
        imageField = possibleImageFields{fieldIdx};
        imgNames = {segData{1}.(imageField)};
        
        % find the name of the other field
        otherField = setdiff(fieldnames(segData{1}), imageField);
        if isempty(otherField)
            error([mfilename ':noseg'], 'Can''t find segmentation data.');
        end
        if numel(otherField) > 1
            error([mfilename ':manyseg'], 'Ambiguous segmentation structure.');
        end
        
        otherField = otherField{1};
        segs = arrayfun(@(s) s.(otherField), segData{1}, 'uniform', false);
        
        if numel(segData) == 2
            segSel = segData{2};
        else
            segSel = [];
        end
    end
end

% set default segSel if it's not provided
initialSegSel = segSel;
if isempty(segSel)
    segSel = ones(numel(segs), 1);
end

% add .jpg to image names that don't have extensions
for i = 1:numel(imgNames)
    [~, ~, ext] = fileparts(imgNames{i});
    if isempty(ext)
        imgNames{i} = [imgNames{i} '.jpg'];
    end
end

% create the main app window
f = figure('Visible', 'off', 'Units', 'normalized', 'Position', settings.figPos);

% make sure segSel is valid
segSel = fixSegChoices(segs, segSel);

% make sure that the field we're trying to use actually exists
settings.field = fixField(settings.field, segs);

% focus on images with segmentations and find their sizes
imageDetails = findImageDetails(segs, settings);

% store common data with the figure
guidata(f, struct(...
    'segs', {segs}, 'segSel', {segSel}, 'settings', settings, ...
    'titleHandles', [], ...
    'imgDirectory', imgDirectory, 'imgNames', {imgNames}, ...
    'imageDetails', {imageDetails},  'displayRow', 1, 'imageAxes', {[]}, ...
    'imagePanels', {[]}, ...
    'positionHandles',  {[]}, 'fieldDropdown', [], ...
    'zoomHandles', [], 'resizeTimer', [], 'controlsHandles', [], ...
    'redrawOnly', [] ...
));

% figure out how to split images by rows
findThumbnailLocations(f);

% create axes for the images
% this also tweaks the vertical positions based on the locations vector
% that was generated by findThumbnailLocations (and added to guidata)
updateAxes(f);

% display controls, including position in image set, buttons to move
% forward and backward, and dropdown menu to select segmentation field
showControls(f);

% make the figure visible
f.Visible = 'on';

% load the images and display them
displayImages(f);

% wait until the fiure is closed
% set a callback for when the figure is about to get closed
set(f, 'closeRequestFcn', @prepareClose, 'resizeFcn', @handleResize);
uiwait(f);

% the 'prepareClose' function didn't allow the figure to get close, which
% allows us to still access guidata, and return the latest version of segSel
data = guidata(f);
segSel = data.segSel;
settings = data.settings;

% save the segmentation if necessary
if fromFile && settings.saveChanges && ~isequaln(segSel, initialSegSel)
    doSave = true;
    if settings.confirmSave
        userChoice = questdlg('There are changes in the segmentation choices. Save to file?', ...
            'Overwrite file?', 'Yes', 'No', 'Yes');
        if ~strcmp(userChoice, 'Yes')
            doSave = false;
        end
    end
    if doSave
        segFileData.(segSelFieldNameFile) = segSel;
        save(segFileName, '-struct', 'segFileData');
    end
end

% save settings for next time the app is used
set(f, 'Units', 'normalized');
settings.figPos = get(f, 'position');
save('segGuiSettings.mat', 'settings');

% now we actually delete the figure
delete(f);

end

function prepareClose(f, ~)
% prepareClose When a figure close is requested, make sure to give the app
% the opportunity to store the output.

uiresume(f);

end

function handleResize(f, ~)
% handleResize Handle a figure resize.

minDelay = 0.5; % minimum delay in seconds before image is redrawn

    function redraw(varargin)
        % need to invalidate control handles!
        sdata = guidata(f);
        
        findThumbnailLocations(f);
        % make sure we're focusing in roughly the same area as before
        idx = sdata.locations{sdata.displayRow}(1).idx;
        
        sdata = guidata(f);
        sdata.displayRow = findRow(f, idx);
        guidata(f, sdata);
        
        updateAxes(f);
        showControls(f);
        displayImages(f);        
    end

data = guidata(f);
if ~isempty(data.resizeTimer)
    stop(data.resizeTimer);
    delete(data.resizeTimer);
end
data.resizeTimer = timer('timerfcn', @redraw, 'startdelay', minDelay);
start(data.resizeTimer);

guidata(f, data);

end

function details = findImageDetails(segs, settings)
% findImageDetails Get image details (like width & height) for those images
% that have segmentations. Return in a structure with a mapping back to
% indices in the original structure.
%
% The field that is used is taken form settings.field.

details = [];

for i = 1:numel(segs)
    crtSeg = segs{i}.(settings.field);
    if ~isempty(crtSeg)
        % fliplr because size(crtSeg, 1) is number of rows, but I want the
        % details.size field to start with the width
        crtDetails = struct('size', fliplr(size(crtSeg)), 'idx', i);
        details = [details ; crtDetails]; %#ok<AGROW>
    end
end

end

function findThumbnailLocations(f)
% findThumbnailLocations Find the horizontal positions for the thumbnails,
% and how to split them by rows.
%   findThumbnailLocations(f) places images horizontally within rows in
%   order to best use the available space on figure f. A field 'locations'
%   is added to the figure's guidata. This field is a cell array of
%   structure arrays; each cell corresponds to one row of images. Each
%   image in a row contains the positions in pixels ('position' field) and
%   the sizes in pixels ('size' field) for the thumbnails. The vertical
%   positions are set according to the minRowSpace option, but may be
%   tweaked before display; in particular, these distances are from the top
%   instead of from the bottom of the image, and they start at 0 instead of
%   topEdge.
%
%   The following settings are used:
%     settings.minRowSpace  -- minimum spacing between rows, in character heights
%     settings.minColSpace  -- minimum spacing between rows, in character widths
%     settings.rowHeight    -- minimum row height in character heights
%     settings.maxAspect    -- horizontal size of thumbnails is capped at maxAspect * row height
%     settings.figEdge      -- figure edge in character widths
%     settings.leftEdge     -- extra edge at figure left (in character widths)
%     settings.rightEdge    -- extra edge at figure right (in character widths)
%     settings.controlsHeight -- space reserved for picture controls (in character heights)
%     settings.titleHeight  -- height reserved for image title (in character heights)

% XXX this function can be quite slow, taking up to two seconds to go
% through 3000 images if only one fits per row

tic;
data = guidata(f);
settings = data.settings;

[charsToPx, sizePx] = getCharsToPx(f);

% figure out how much of the figure is taken by edges on each side
figEdgePx = round(settings.figEdge * charsToPx(1));
extraLeftEdgePx = round(settings.leftEdge * charsToPx(1));
extraRightEdgePx = round(settings.rightEdge * charsToPx(1));

leftEdgePx = figEdgePx + extraLeftEdgePx;
rightEdgePx = figEdgePx + extraRightEdgePx;

% figure out row height and distance between rows
rowHeightPx = round(settings.rowHeight*charsToPx(2));
minRowSpacePx = settings.minRowSpace*charsToPx(2);
controlsHeightPx = round(settings.controlsHeight*charsToPx(2));
titleHeightPx = settings.titleHeight*charsToPx(2);

minTotalRowHeightPx = rowHeightPx + titleHeightPx + controlsHeightPx + minRowSpacePx;

% figure out minimum horizontal spacing
minColSpacePx = round(settings.minColSpace * charsToPx(1));

% figure out thumbnail sizes
imageDetails = data.imageDetails;
imageWidths = arrayfun(@(d) d.size(1), imageDetails);
imageHeights = arrayfun(@(d) d.size(2), imageDetails);
thumbHeights = repmat(rowHeightPx, size(imageWidths));
thumbWidths = imageWidths*rowHeightPx ./imageHeights;
% find images that are too wide and scale them down
maxWidthPx = round(rowHeightPx*settings.maxAspect);
maskTooWide = (thumbWidths > maxWidthPx);
thumbHeights(maskTooWide) = round(imageHeights(maskTooWide)*maxWidthPx ...
    ./ imageWidths(maskTooWide));
thumbWidths(maskTooWide) = maxWidthPx;

% estimate number of rows needed
cumulativeWidths = cumsum(thumbWidths + minColSpacePx);
estimatedRows = ceil(cumulativeWidths(end) / ...
    (sizePx(1) - leftEdgePx - rightEdgePx));

% start placing the figures
locations = cell(estimatedRows, 1);
imgIdx = 1;
rowIdx = 1;
nImg = numel(imageDetails);
rowThresholdPx = sizePx(1) - leftEdgePx - rightEdgePx + minColSpacePx;
flatten = @(x) x(:);
while imgIdx <= nImg
    % find the last image that fits on this row
    if imgIdx > 1
        startCum = cumulativeWidths(imgIdx - 1);
    else
        startCum = 0;
    end
    nImagesPlusOne = find(...
        cumulativeWidths(imgIdx:end) - startCum > rowThresholdPx, 1);
    if isempty(nImagesPlusOne)
        lastImgIdx = nImg;
    else
        if nImagesPlusOne > 1
            lastImgIdx = imgIdx + nImagesPlusOne - 2;
        else
            % make sure there's at least one image per row
            lastImgIdx = imgIdx;
        end
    end
    nImages = lastImgIdx - imgIdx + 1;
    rowImagesWidth = sum(thumbWidths(imgIdx:lastImgIdx));
    if nImages > 1
        imageSpacingPx = (sizePx(1) - leftEdgePx - rightEdgePx - rowImagesWidth) ...
            / (nImages - 1);
    else
        imageSpacingPx = 0;
    end
    
    xLocations = [0 ; flatten(cumsum(thumbWidths(imgIdx:lastImgIdx)))];
    locationsRow = struct(...
        'size', arrayfun(@(i) [thumbWidths(i), thumbHeights(i)], imgIdx:lastImgIdx, 'uniform', false), ...
        'idx', num2cell(imgIdx:lastImgIdx), ...
        'row', repmat(rowIdx, nImages, 1), ...
        'position', arrayfun(@(i) ...
            [xLocations(i - imgIdx + 1) + (i-imgIdx)*imageSpacingPx + leftEdgePx, ...
                (rowIdx - 1)*minTotalRowHeightPx], ...
            imgIdx:lastImgIdx, 'uniform', false) ...
    );

    locations{rowIdx} = locationsRow;
    
    imgIdx = lastImgIdx + 1;
    rowIdx = rowIdx + 1;
end
nRows = rowIdx - 1;
if numel(locations) > nRows
    locations = locations(1:nRows);
end

data.locations = locations;
guidata(f, data);
%disp(['Find thumbnail locations took ' num2str(toc, '%.3f') ' seconds.']);

end

function s = setIfMissing(s, varargin)
% setIfMissing Set structure fields only if they don't yet exist or if they
% are empty.

keys = varargin(1:2:end);
values = varargin(2:2:end);

for i = 1:length(keys)
    if ~isfield(s, keys{i}) || isempty(s.(keys{i}))
        s.(keys{i}) = values{i};
    end
end

end

function updateAxes(f)
% updateAxes Create axes at the locations (in pixels) given in the figure's
% guidata.
%   updateAxes(f) creates axes for images identified from the 'locations'
%   field in the figure's guidata. This starts at row 'displayRow' from the
%   guidata, and uses the size of the figure 'f' and the settings (from
%   guidata) to decide how many rows of images to display. This is saved in
%   guidata as 'nRows'. The function also either deletes any old axes from
%   the guidata field 'imageAxes', and replaces them with the new ones; or
%   reuses the old axes.
%
%   Settings used
%     settings.minRowSpace  -- minimum spacing between rows, in character heights
%     settings.minColSpace  -- minimum spacing between rows, in character widths
%     settings.rowHeight    -- minimum row height in character heights
%     settings.topEdge      -- extra edge at top of figure (in character heights)
%     settings.bottomEdge   -- extra edge at bottom of figure (in character heights)
%     settings.controlsHeight -- space reserved for picture controls (in character heights)
%     settings.titleHeight  -- height reserved for image title (in character heights)

data = guidata(f);
ax = data.imageAxes;
panels = data.imagePanels;
locations = data.locations;
row0 = data.displayRow;
settings = data.settings;

[charsToPx, sizePx] = getCharsToPx(f);

% figure out how much of the figure is taken by edges
figEdgePx = round(settings.figEdge * charsToPx(1));

extraTopEdgePx = round(settings.topEdge * charsToPx(2));
extraBottomEdgePx = round(settings.bottomEdge * charsToPx(2));

topEdgePx = figEdgePx + extraTopEdgePx;
bottomEdgePx = figEdgePx + extraBottomEdgePx;

% figure out how much vertical space each thumbnail takes
rowHeightPx = round(settings.rowHeight*charsToPx(2));
minRowSpacePx = settings.minRowSpace*charsToPx(2);
controlsHeightPx = round(settings.controlsHeight*charsToPx(2));
titleHeightPx = settings.titleHeight*charsToPx(2);

fixedHeightPx = rowHeightPx + titleHeightPx + controlsHeightPx;
minTotalRowHeightPx = fixedHeightPx + minRowSpacePx;

% useful figure height
usefulHeightPx = sizePx(2) - topEdgePx - bottomEdgePx;
nRows = floor((usefulHeightPx + minRowSpacePx) / minTotalRowHeightPx);
if nRows > 1
    actualRowSkipPx = (usefulHeightPx - fixedHeightPx) / (nRows - 1);
else
    actualRowSkipPx = 0;
    nRows = 1;
end

crtY = sizePx(2) - topEdgePx - titleHeightPx;
axIdx = 1;
for j = 1:nRows
    row = row0 + j - 1;
    if row > numel(locations)
        break;
    end
    crtLocations = locations{row};
    for i = 1:numel(crtLocations)
        % if the thumbnail is smaller than usual, center it
        yPos = crtY - 0.5*(crtLocations(i).size(2) + rowHeightPx);
        
        args0 = {'units', 'pixels', 'position', ...
            [0 0 round(crtLocations(i).size)]};
        args = {'units', 'pixels', 'position', ...
            round([crtLocations(i).position(1) yPos ...
             crtLocations(i).size])};
        if axIdx <= numel(ax)
            % reuse existing axes
            set(panels(axIdx), args{:});
            set(ax(axIdx), args0{:});
        else
            % or create new ones if needed
            panels(axIdx) = uipanel('bordertype', 'none', args{:});
            ax(axIdx) = axes(args0{:},'parent', panels(axIdx));
        end
        
        set(ax, 'xtick', [], 'ytick', []);
        
        axIdx = axIdx + 1;
    end
    
    crtY = crtY - actualRowSkipPx;
end
nAxes = axIdx - 1;
% delete axes that we no longer need
if numel(ax) > nAxes
    delete(ax(nAxes+1:end));
    ax = ax(1:nAxes);
    
    delete(panels(nAxes+1:end));
    panels = panels(1:nAxes);
end

data.imageAxes = ax;
data.imagePanels = panels;
data.nRows = nRows;
guidata(f, data);

end

function segSel = fixSegChoices(segs, segSel)
% fixSegChoices Make sure segmentation choices are valid.

for i = 1:numel(segs)
    crtSeg = segs{i};
    crtSegChoice = segSel(i);
    if isfinite(crtSegChoice)
        if crtSegChoice < 0 || crtSegChoice > numel(crtSeg)
            crtSegChoice = 1;
            segSel(i) = crtSegChoice;
        end
    end
end

end

function displayImages(f)
% displayImages Display the images in the axes.
%   The handles for the images are stored in the 'imageHandles' field of
%   the figure's guidata.

data = guidata(f);
ax = data.imageAxes;
panels = data.imagePanels;
locations = data.locations;
row0 = data.displayRow;
details = data.imageDetails;
segs = data.segs;
segSel = data.segSel;
directory = data.imgDirectory;
names = data.imgNames;
showSeg = data.settings.showSeg;
%settings = data.settings;

segField = data.settings.field;
titleHandles = data.titleHandles;
controlsHandles = data.controlsHandles;
redrawOnly = data.redrawOnly;

charsToPx = getCharsToPx(f);

n = length(ax);

h = [];

row = row0;
withinRow = 1;
for i = 1:n
    % find the segmentation
    crtLocation = locations{row}(withinRow);
    
    segIdx = details(crtLocation.idx).idx;
    
    if isempty(redrawOnly) || ismember(segIdx, redrawOnly)
        crtSegAll = segs{segIdx};
        crtSegChoice = segSel(segIdx);
        crtSeg = crtSegAll;
        if isfinite(crtSegChoice)
            crtSeg = crtSeg(crtSegChoice);
        elseif ~isempty(crtSeg)
            crtSeg = crtSeg(1);
        end
        
        % find the image
        crtName = fullfile(directory, names{segIdx});
        crtImage = loadLUMImage(crtName);
        
        axPos = get(panels(i), 'position');
        if showSeg
            % show segmentation, unless the user chose not to
            
            % check that the image and the segmentation match
            if ~isempty(segField) && ~all(size(crtSeg.(segField)) == size(crtImage))
                warning([mfilename ':segsizemis'], 'Segmentation size not matching image size.');
            end
            
            % sanity check that the image and the axes in which it is to be drawn match
            aspectImage = size(crtImage, 2) / size(crtImage, 1);
            aspectLoc = axPos(3) / axPos(4);
            relErrImage = sum(0.5 ./ size(crtImage));
            relErrLoc = sum(0.5 ./ axPos(3:4));
            if abs(aspectLoc - aspectImage) / mean([aspectLoc, aspectImage]) > relErrLoc + relErrImage
                warning([mfilename ':locsizemis'], 'Location aspect ratio not matching image.');
            end
        end
        
        % display the image
        [crtH, crtHasSeg, crtIgnoreSeg] = displaySingleImage(...
            ax(i), crtImage, crtSeg, crtSegChoice, segField, showSeg);
        
        if ~showSeg
            crtHasSeg = false;
        end
        
        % set a callback to detect double click
%        set(ax(i), 'buttondownfcn', {@thumbClickHandler, f, crtLocation});
        set(crtH, 'buttondownfcn', {@thumbClickHandler, f, crtLocation});
                
        % display image title
        titleArgs = {'position', round([...
            axPos(1) axPos(2) + axPos(4) axPos(3) charsToPx(2)]), ...
            'string', names{segIdx}, 'tooltipstring', names{segIdx}};
        if i <= numel(titleHandles)
            % reuse existing handles
            set(titleHandles(i), titleArgs{:});
        else
            % or create a new one
            titleHandles(i) = uicontrol('style', 'text', titleArgs{:});
        end
        titleExtent = get(titleHandles(i), 'extent');
        if titleExtent(3) > axPos(3)
            % the title didn't fit, display a short version with '...' at the end
            set(titleHandles(i), 'string', [names{segIdx}(1:3) '...']);
        end
        
        % display the controls for each image
        ctrlIdx = 1 + (i-1)*4;
        argsPrev = {'position', round([...
            1 1 charsToPx(2) charsToPx(2)]), ...
            'string', '<', 'tooltipstring', 'Previous segmentation'};
        argsNoSeg = {'position', round([...
            1+charsToPx(2) 1 charsToPx(2) charsToPx(2)]), ...
            'string', 'x', 'tooltipstring', 'Ignore this image'};
        argsNext = {'position', round([...
            1+2*charsToPx(2) 1 charsToPx(2) charsToPx(2)]), ...
            'string', '>', 'tooltipstring', 'Next segmentation'};
        argsPanel = {'units', 'pixels', 'position', round([...
            axPos(1) axPos(2)-charsToPx(2)-2 3*charsToPx(2)+2 charsToPx(2)])};
        if ctrlIdx < numel(controlsHandles)
            % reuse existing handles
            set(controlsHandles(ctrlIdx+3), argsPanel{:});
            set(controlsHandles(ctrlIdx), argsPrev{:});
            set(controlsHandles(ctrlIdx+1), argsNoSeg{:});
            set(controlsHandles(ctrlIdx+2), argsNext{:});
        else
            % create new panel
            controlsHandles(ctrlIdx+3) = uipanel('bordertype', 'none', ...
                argsPanel{:});
            
            % create new buttons
            controlsHandles(ctrlIdx) = uicontrol('parent', controlsHandles(ctrlIdx+3), ...
                'style', 'pushbutton', argsPrev{:});
            controlsHandles(ctrlIdx+1) = uicontrol('parent', controlsHandles(ctrlIdx+3), ...
                'style', 'pushbutton', argsNoSeg{:});
            controlsHandles(ctrlIdx+2) = uicontrol('parent', controlsHandles(ctrlIdx+3), ...
                'style', 'pushbutton', argsNext{:});
        end
        
        % enable or disable arrow buttons as needed
        set(controlsHandles(ctrlIdx+1), 'enable', 'on');
        if ~crtHasSeg || crtIgnoreSeg
            set(controlsHandles([ctrlIdx, ctrlIdx+2]), 'enable', 'off');
        else
            set(controlsHandles([ctrlIdx, ctrlIdx+2]), 'enable', 'on');
        end
        if crtHasSeg && ~crtIgnoreSeg
            if crtSegChoice == 1
                set(controlsHandles(ctrlIdx), 'enable', 'off');
            end
            if crtSegChoice == numel(crtSegAll)
                set(controlsHandles(ctrlIdx+2), 'enable', 'off');
            end
        end
        
        set(controlsHandles(ctrlIdx+1), 'callback', {@toggleIgnoreImage, f, segIdx});
        
        set(controlsHandles(ctrlIdx), 'callback', {@shiftSegSel, f, segIdx, -1});
        set(controlsHandles(ctrlIdx+2), 'callback', {@shiftSegSel, f, segIdx, 1});
        
        % store the handles for everything, so they can be easily either deleted
        % or reused when an update is necessary
        h = [h ; crtH]; %#ok<AGROW>
    end
    
    % move along within the row, or on to the next row
    withinRow = withinRow + 1;
    if withinRow > numel(locations{row})
        % redraw once per row
        %drawnow('limitrate', 'nocallbacks');
        row = row + 1;
        withinRow = 1;
        if row > numel(locations)
            break;
        end
    end
end

% cull un-needed handles
if numel(titleHandles) > n
    delete(titleHandles(n+1:end));
    titleHandles = titleHandles(1:n);
end
if numel(controlsHandles) > 4*n
    delete(controlsHandles(4*n+1:end));
    controlsHandles = controlsHandles(1:4*n);
end

% save the handles
data.titleHandles = titleHandles;
if isempty(redrawOnly)
    data.imageHandles = h;
end
data.controlsHandles = controlsHandles;
data.redrawOnly = [];
guidata(f, data);

% redraw
drawnow('limitrate');

end

function showControls(f)
% showControls Display controls to move around in the image set and to
% select segmentation field.
%   showControls(f) displays
%       - a label showing which images are currently visible
%       - buttons to move backward and foward through the images
%       - a dropdown menu to choose the segmentation field to use
%       - buttons to increase or decrease the size of the thumbnails
%   The handles for the position label and buttons are stored (in the
%   order label, button prev, button next) in the vector 'positionHandles'
%   stored in the figure's guidata. The handles for the dropdown and its
%   label are stored in the 'fieldDropdown' field (in the order
%   [label, dropdown]). The handles for the '+' and '-' zoom buttons are
%   stored in the field 'zoomHandles'.
%
%   The function uses the fields 'locations', 'displayRow', 'nRows', and
%   'settings' from the guidata. If the handles are already set to non-NaN
%   non-empty values, the controls are simply updated.

data = guidata(f);
h = data.positionHandles;
if isempty(h)
    h = [nan, nan, nan, nan];
end
locations = data.locations;
row0 = data.displayRow;
row1 = row0 + data.nRows - 1;
settings = data.settings;

[charsToPx, sizePx] = getCharsToPx(f);

% figure out how much of the figure is taken by edges on each side
figEdgePx = round(settings.figEdge * charsToPx(1));
extraRightEdgePx = round(settings.rightEdge * charsToPx(1));
rightEdgePx = figEdgePx + extraRightEdgePx;

% display the label showing the range of thumbnails we're looking at
maxPositionLength = 20;
labelWidthPx = maxPositionLength*charsToPx(1);
if isnan(h(1))    
    h(1) = uicontrol('style', 'text', 'horizontalalignment', 'right', ...
        'tooltipstring', 'Double-click to select index', 'enable', 'inactive');
end

nTotal = locations{end}(end).idx;
if row1 > numel(locations)
    row1 = numel(locations);
end
set(h(1), 'string', [int2str(locations{row0}(1).idx) ' - ' ...
    int2str(locations{row1}(end).idx) ' of ' int2str(nTotal)], ...
    'position', round([...
        sizePx(1) - rightEdgePx - labelWidthPx, ...
        sizePx(2) - figEdgePx - charsToPx(2), ...
        labelWidthPx, charsToPx(2)]));
set(h(1), 'buttondownfcn', {@gotoIndexClickHandler, f, locations{row0}(1).idx});

% display the 'prev' and 'next' buttons
maxButtonLength = 8;
buttonLengthPx = maxButtonLength*charsToPx(1);
interButtonDistance = charsToPx(1);
totalButtonLength = 2*buttonLengthPx + interButtonDistance;

% make the panel containing the 'prev' and 'next' buttons
panelX = sizePx(1) - rightEdgePx - totalButtonLength;
panelY = figEdgePx/2;
panelWidth = totalButtonLength;
panelHeight = 2*charsToPx(2);
if isnan(h(4))
    h(4) = uipanel('bordertype', 'none');
end
set(h(4), 'units', 'pixels', ...
    'position', round([panelX panelY panelWidth panelHeight]));

if isnan(h(2))
    h(2) = uicontrol('parent', h(4), 'style', 'pushbutton', 'string', 'Prev', ...
        'tooltipstring', 'Show previous page');
end
set(h(2), 'units', 'pixels', ...
    'position', round([0, 0, buttonLengthPx, panelHeight]));

if isnan(h(3))
    h(3) = uicontrol('parent', h(4), 'style', 'pushbutton', 'string', 'Next', ...
        'tooltipstring', 'Show next page');
end
set(h(3), 'units', 'pixels', ...
    'position', round([...
            buttonLengthPx + interButtonDistance, 0, buttonLengthPx, panelHeight]));

% update the button's visibility and callbacks
if row0 <= 1
    set(h(2), 'enable', 'off');
else
    set(h(2), 'enable', 'on');
    set(h(2), 'callback', {@jumpToRow, f, max(1, row0-data.nRows)});
end

if row1 >= numel(locations)
    set(h(3), 'enable', 'off');
else
    set(h(3), 'enable', 'on');
    set(h(3), 'callback', {@jumpToRow, f, min(numel(locations), row0+data.nRows)});
end

% display or update the dropdown
fieldLabelText = 'field:';
fieldLabelWidth = charsToPx(1)*5;
fieldToggleWidth = 1.3*charsToPx(2);
fieldLabelX = sizePx(1)/2 - fieldLabelWidth;
fieldLabelY = sizePx(2) - figEdgePx - charsToPx(2);
fieldDropdownWidth = charsToPx(1)*16;

if isempty(data.fieldDropdown)
    fieldPanel = uipanel('bordertype', 'none');
    fieldLabelHandle = uicontrol('parent', fieldPanel, 'style', 'text', ...
        'string', fieldLabelText);
    
    fieldDropdownHandle = uicontrol('parent', fieldPanel, 'style', 'popup');
    fieldToggleHandle = uicontrol('parent', fieldPanel, 'style', 'checkbox');
    
    data.fieldDropdown = [fieldLabelHandle fieldDropdownHandle fieldPanel fieldToggleHandle];
end
set(data.fieldDropdown(3), 'units', 'pixels', ...
    'position', round([fieldLabelX fieldLabelY-charsToPx(2) ...
        fieldLabelWidth + fieldDropdownWidth + fieldToggleWidth charsToPx(2)*4]));
set(data.fieldDropdown(1), 'position', round([fieldToggleWidth, charsToPx(2), ...
            fieldLabelWidth, charsToPx(2) ...
        ]));
set(data.fieldDropdown(2), 'position', round([...
            fieldToggleWidth + fieldLabelWidth, 0.2*charsToPx(2), ...
            fieldDropdownWidth, charsToPx(2)*2 ...
        ]));
set(data.fieldDropdown(4), 'position', round([...
    0, charsToPx(2), fieldToggleWidth, 1.1*charsToPx(2)]));
checkMinMax = [get(data.fieldDropdown(4), 'min') get(data.fieldDropdown(4), 'max')];
set(data.fieldDropdown(4), 'value', checkMinMax(data.settings.showSeg + 1), ...
    'callback', {@changeSegVisibility, f});
if ~data.settings.showSeg
    set(data.fieldDropdown(1:2), 'enable', 'off');
else
    set(data.fieldDropdown(1:2), 'enable', 'on');
end

% figure out what fields we have
segs = data.segs;
segNames = fieldnames(segs{1});
% filter out fields that don't look like segmentation masks
isGood = @(m) ismatrix(m) && ~isvector(m) && isnumeric(m);
mask = cellfun(@(name) isGood(segs{1}(1).(name)), segNames);
segNames = segNames(mask);

% add a 'none' field -- don't display any segmentation
% segNames = [{'none'} ; segNames(:)];

crtSegIdx = find(strcmp(segNames, data.settings.field));
set(data.fieldDropdown(2), 'string', segNames, 'value', crtSegIdx, ...
    'callback', {@changeField, f});

% display the zoom buttons
if isempty(data.zoomHandles)
    minusHandle = uicontrol('style', 'pushbutton', 'string', '-', ...
        'tooltipstring', 'Decrease thumbnail size');
    plusHandle = uicontrol('style', 'pushbutton', 'string', '+', ...
        'tooltipstring', 'Increase thumbnail size');
    
    data.zoomHandles = [minusHandle plusHandle];
end
set(data.zoomHandles(1), 'position', round([figEdgePx, sizePx(2) - figEdgePx - charsToPx(2), ...
                     charsToPx(2), charsToPx(2)]));
set(data.zoomHandles(2), 'position', round([figEdgePx + charsToPx(2), sizePx(2) - figEdgePx - charsToPx(2), ...
                     charsToPx(2), charsToPx(2)]));

minusFactor = 3/4;
plusFactor = 4/3;
if validateZoomFactor(f, minusFactor)
    set(data.zoomHandles(1), 'enable', 'on');
    set(data.zoomHandles(1), 'callback', {@zoomThumbnails, f, minusFactor});
else
    set(data.zoomHandles(1), 'enable', 'off');
end

if validateZoomFactor(f, plusFactor)
    set(data.zoomHandles(2), 'enable', 'on');
    set(data.zoomHandles(2), 'callback', {@zoomThumbnails, f, plusFactor});
else
    set(data.zoomHandles(2), 'enable', 'off');
end

% update guidata
data.positionHandles = h;
guidata(f, data);

end

function jumpToRow(~, ~, f, row)
% jumpToRow Display images starting at the given row.

data = guidata(f);
data.displayRow = row;

guidata(f, data);

updateAxes(f);
showControls(f);
displayImages(f);

end

function changeField(dropdown, ~, f)
% changeField Change the segmentation field that is currently in use.

idx = dropdown.Value;
segNames = dropdown.String;

newField = segNames{idx};

data = guidata(f);
data.settings.field = newField;
guidata(f, data);

%updateAxes(f);
showControls(f);
displayImages(f);

end

function changeSegVisibility(checkbox, ~, f)
% changeSegVisibility Change visibility of segmentations.

data = guidata(f);
data.settings.showSeg = (get(checkbox, 'value') == get(checkbox, 'max'));
guidata(f, data);

showControls(f);
drawnow('limitrate');
displayImages(f);

end

function field = fixField(field, segs)
% fixField Check whether given field name is valid; if not, choose first
% segmentation-like field from the 'segs' structure.

if ~isfield(segs{1}, field)
    oldField = field;
    
    segNames = fieldnames(segs{1});
    % filter out fields that don't look like segmentation masks
    isGood = @(m) ismatrix(m) && ~isvector(m) && isnumeric(m);
    mask = cellfun(@(name) isGood(segs{1}(1).(name)), segNames);
    segNames = segNames(mask);
    
    if isempty(segNames)
        field = '';
        warning('Segmentation structure does not have any valid fields.');
    else
        field = segNames{1};
        warning(['Segmentation structure does not have field ' oldField '; using ' field ' instead.']);
    end
end

end

function zoomThumbnails(~, ~, f, factor)
% zoomThumbnails Rescale the size of the thumbnails by the given factor.

data = guidata(f);

[valid, newRowHeight] = validateZoomFactor(f, factor);
if ~valid
    return;
end

data.settings.rowHeight = newRowHeight;

guidata(f, data);

findThumbnailLocations(f);
% make sure we're focusing in roughly the same area as before
idx = data.locations{data.displayRow}(1).idx;

data = guidata(f);
data.displayRow = findRow(f, idx);
guidata(f, data);

updateAxes(f);
showControls(f);
displayImages(f);

end

function [valid, newRowHeight] = validateZoomFactor(f, factor)
% validateZoomFactor Check whether the given zoom factor leads to a
% reasonable rowHeight.

data = guidata(f);
newRowHeight = data.settings.rowHeight * factor;

valid = false;
if newRowHeight < 3
    % too small
    newRowHeight = data.settings.rowHeight;
    return;
end

[~, ~, sizeChars] = getCharsToPx(f);
if newRowHeight >= sizeChars(2)
    % too big
    newRowHeight = data.settings.rowHeight;
    return;
end

valid = true;

end

function [charsToPx, sizePx, sizeChars] = getCharsToPx(f)
% getCharsToPx Get the conversion from character width/height to pixels.

% get size of window in pixels and in characters
oldUnits = get(f, 'units');

set(f, 'units', 'pixels');
positionPx = get(f, 'position');
sizePx = positionPx(3:4);

set(f, 'units', 'characters');
positionChars = get(f, 'position');
sizeChars = positionChars(3:4);
% get conversion from character widths/heights to pixels
charsToPx = sizePx ./ sizeChars;

set(f, 'units', oldUnits);

end

function row = findRow(f, idx)
% findRow Find the row in which the given thumbnail is found.

data = guidata(f);

idxs = cellfun(@(locRow) arrayfun(@(l) l.idx, locRow), data.locations, ...
    'uniform', false);
mask = cellfun(@(ix) ismember(idx, ix), idxs);
row = find(mask, 1);

end

function toggleIgnoreImage(~, ~, f, idx)
% toggleIgnoreImage Change the selection for the given image to either NaN
% (representing ignoring), or from NaN to 1 if it was already ignored.

data = guidata(f);
if isnan(data.segSel(idx))
    data.segSel(idx) = 1;
else
    data.segSel(idx) = nan;
end

data.redrawOnly = idx;

guidata(f, data);

displayImages(f);

end

function shiftSegSel(~, ~, f, idx, shift)
% shiftSegSel Shift the segmentation selection for the given image. This
% assumes that the given image is not set to NaN, and also that the shift
% does not lead to an invalid segmentation index.

data = guidata(f);
data.segSel(idx) = data.segSel(idx) + shift;

data.redrawOnly = idx;
guidata(f, data);

displayImages(f);

end

function [crtH, crtHasSeg, crtIgnoreSeg] = displaySingleImage(ax, crtImage, ...
    crtSeg, crtSegChoice, segField, showSeg)
% displaySingleImage Display the image in the given axes, overlaying the
% appropriate segmentation (or a cross if crtSegChoice is NaN). Set
% 'showSeg' to false to not display the segmentation.

% clear the axes
delete(get(ax, 'children'));

% turn the image into grayscale RGB, normalizing intensities
normalizedImage = crtImage / max(crtImage(:));

% overlay the segmentation, if available
overlayPlanes = {normalizedImage, normalizedImage, normalizedImage};
crtHasSeg = ~isempty(segField) && ~strcmp(segField, 'none') && ~isempty(crtSeg);
crtIgnoreSeg = isnan(crtSegChoice);
%    if isfinite(crtSegChoice) && ~isempty(crtSeg) && ~isempty(segField)
if showSeg && crtHasSeg % && ~crtIgnoreSeg
    crtSegImage = crtSeg.(segField);
    segmentValues = setdiff(unique(crtSegImage(:)), 0);
    nSegments = numel(segmentValues);
    segColors = parula(nSegments);
    
    for j = 1:nSegments
        crtColor = segColors(j, :);
        for k = 1:3
            mask = crtSegImage == segmentValues(j);
            overlayPlanes{k}(mask) = 0.25*overlayPlanes{k}(mask) + 0.75*crtColor(k);
        end
    end
end
%    if isnan(crtSegChoice)
if crtIgnoreSeg
    % this image is to be ignored
    for k = 1:3
        overlayPlanes{k} = 0.5*overlayPlanes{k};
    end
end
overlay = zeros([size(normalizedImage) 3]);
for k = 1:3
    overlay(:, :, k) = overlayPlanes{k};
end

% draw the overlay
crtH = image(uint8(255*overlay), 'parent', ax);
colormap('gray');
%    if isnan(crtSegChoice)
if crtIgnoreSeg
    % this image is to be ignored -- draw an X over it
    lineHandle1 = line([0 size(normalizedImage, 2)], ...
        [0 size(normalizedImage, 1)], 'parent', ax);
    lineHandle2 = line([size(normalizedImage, 2) 0], ...
        [0 size(normalizedImage, 1)], 'parent', ax);
    set([lineHandle1, lineHandle2], 'linewidth', 2, 'color', 'r');
end
% get rid of ticks
set(ax, 'xtick', [], 'ytick', []);

end

function gotoIndexClickHandler(~, ~, f, crtIdx)
% gotoIndexClickHandler Handle double click on label showing position
% within image set.

persistent clickCount;

if isempty(clickCount)
    clickCount = 1;
    pause(0.4);
    if clickCount == 1
        clickCount = [];
    end
else
    newIdxStr = inputdlg({'Go to image index:'}, 'Go to image', 1, {int2str(crtIdx)});
    try
        newIdx = str2double(newIdxStr);
        if newIdx == round(newIdx)
            newRow = findRow(f, newIdx);
            jumpToRow(0, 0, f, newRow);
        end
    catch
    end
    clickCount = [];
end

end

function thumbClickHandler(~, ~, f, loc)
% thumbClickHandler Handle double click on thumbnail.

persistent clickCount;

if isempty(clickCount)
    clickCount = 1;
    pause(0.4);
    if clickCount == 1
        clickCount = [];
    end
else
    data = guidata(f);
    data = singleImageDialog(data, loc.idx);
    data.redrawOnly = loc.idx;
    guidata(f, data);
    displayImages(f);
    clickCount = [];
end

end

function data = singleImageDialog(data, locIdx)
% singleImageDialog Process a single image.

% create a modal dialog box
segIdx = data.imageDetails(locIdx).idx;
d = dialog('units', 'normalized', 'position', data.settings.singleImagePos, ...
    'name', data.imgNames{segIdx}, 'resize', 'on');

% get the segmentation info and draw the image
crtSegAll = data.segs{segIdx};
crtSegChoice = data.segSel(segIdx);
crtSeg = crtSegAll;
if isfinite(crtSegChoice)
    crtSeg = crtSeg(crtSegChoice);
elseif ~isempty(crtSeg)
    crtSeg = crtSeg(1);
end

% find the image
crtName = fullfile(data.imgDirectory, data.imgNames{segIdx});
crtImage = loadLUMImage(crtName);

% figure out scaling
[charsToPx, sizePx] = getCharsToPx(d);

figEdgePx = round(data.settings.figEdge * charsToPx(1));
factor = min((sizePx - figEdgePx) ./ fliplr(size(crtImage)));

axesSize = round(fliplr(size(crtImage))*factor);

% draw the image
ax = axes('units', 'pixels', 'position', [0.5*(sizePx(1:2) - axesSize(1:2)) ...
    axesSize]);
[hLeftArrow, hRightArrow] = redrawImage;
set(d, 'windowbuttonmotionfcn', @mouseMove, 'windowbuttondownfcn', @mouseDown, ...
    'resizefcn', @dialogResize);
drawnow;

% give us a chance to find out the final position of the dialog before
% deleting it
set(d, 'closerequestfcn', @prepareCloseDlg);

% wait for the dialog box to end
lastDialogSize = [];
uiwait(d);
drawnow;
pause(0.1);

data.segSel(segIdx) = crtSegChoice;
data.settings.singleImagePos = lastDialogSize;

    function prepareCloseDlg(dd, ~)
        set(dd, 'units', 'normalized');
        lastDialogSize = get(dd, 'position');
        delete(dd);
    end

    function h = drawArrow(tipPos, height, orientation, thickness, color, alpha)
        % Draw an arrow in an axis object.
        %   tipPos is a the position of the tip of the arrow in normalized
        %       coordinates.
        %   height is the total height of the arrow in normalized
        %       coordinates.
        %   orientation is 1 for a right-pointing arrow, -1 for a
        %       left-pointing arrow.
        %   thickness gives the width of the line used to draw the arrow.
        %   color gives the color of the arrow.
        
        xRange = xlim;
        yRange = ylim;
        
        tipPosScaled = [xRange(1) yRange(1)] + tipPos .* [diff(xRange) diff(yRange)];
        
        heightScaled = height*diff(yRange);
        widthScaled = -orientation*heightScaled*0.3;
        thicknessScaled = thickness/500*diff(xRange);
        
        edge1x = [tipPosScaled(1) + widthScaled tipPosScaled(1) tipPosScaled(1) + widthScaled];
        edge1y = [tipPosScaled(2) + heightScaled tipPosScaled(2) tipPosScaled(2) - heightScaled];
        h = patch([edge1x(:) ; flipud(edge1x(:)) - orientation*thicknessScaled], ...
            [edge1y(:) ; flipud(edge1y(:))], color);
        
        set(h, 'edgecolor', 'w', 'facealpha', alpha, 'edgealpha', alpha);
    end

    function b = isIn(c, boundsX, boundsY)
        b = (c(1) >= boundsX(1) && c(1) <= boundsX(2) && c(2) >= boundsY(1) && c(2) <= boundsY(2));
    end

    function b = hitsArrow(c, arrow)
        limX = [min(get(arrow, 'xdata')) max(get(arrow, 'xdata'))];
        limY = [min(get(arrow, 'ydata')) max(get(arrow, 'ydata'))];
        b = isIn(c, limX, limY);
    end

    function mouseMove(~, ~)        
        c = get(ax, 'currentpoint');
        c = c(1, :);
        
        xRange = xlim(ax);
        yRange = ylim(ax);
        
        if data.settings.showSeg && ~isnan(crtSegChoice) && isIn(c, xRange, yRange)
            leftVisible = (crtSegChoice > 1);
            rightVisible = (crtSegChoice < numel(crtSegAll));
            
            if leftVisible
                set(hLeftArrow, 'visible', 'on');
                leftAlpha = 0.5 + 0.5*hitsArrow(c, hLeftArrow);
                set(hLeftArrow, 'facealpha', leftAlpha, 'edgealpha', leftAlpha);
            else
                set(hLeftArrow(:), 'visible', 'off');
            end
            
            if rightVisible
                set(hRightArrow, 'visible', 'on');
                rightAlpha = 0.5 + 0.5*hitsArrow(c, hRightArrow);
                set(hRightArrow, 'facealpha', rightAlpha, 'edgealpha', rightAlpha);
            end
        else
            set([hLeftArrow(:) ; hRightArrow], 'visible', 'off');
        end
    end

    function mouseDown(~, ~)
        c = get(ax, 'currentpoint');
        c = c(1, :);
        
        xRange = xlim(ax);
        yRange = ylim(ax);

        if isIn(c, xRange, yRange)
            if data.settings.showSeg
                if ~isnan(crtSegChoice)
                    leftVisible = strcmp(get(hLeftArrow, 'visible'), 'on');
                    rightVisible = strcmp(get(hRightArrow, 'visible'), 'on');
                    
                    if leftVisible && hitsArrow(c, hLeftArrow)
                        shiftSegmentation(-1);
                        return;
                    end
                    if rightVisible && hitsArrow(c, hRightArrow)
                        shiftSegmentation(1);
                        return;
                    end
                    
                    % if we're not clicking on the arrows, we're clicking
                    % directly on the image --> toggle whether to ignore image
                    crtSegChoice = nan;
                else
                    crtSegChoice = 1;
                end
            else
                if isnan(crtSegChoice)
                    crtSegChoice = 1;
                else
                    crtSegChoice = nan;
                end
            end
            [hLeftArrow, hRightArrow] = redrawImage;
        end
    end

    function shiftSegmentation(shift)
        crtSegChoice = crtSegChoice + shift;
        crtSeg = crtSegAll;
        if isfinite(crtSegChoice)
            crtSeg = crtSeg(crtSegChoice);
        elseif ~isempty(crtSeg)
            crtSeg = crtSeg(1);
        end
        
        [hLeftArrow, hRightArrow] = redrawImage;
    end

    function [hl, hr] = redrawImage
        displaySingleImage(ax, crtImage, crtSeg, crtSegChoice, data.settings.field, ...
            data.settings.showSeg);
        hl = drawArrow([0.05 0.5], 0.3, -1, 5, [0.8 0.4 0.4], 0.5);
        hr = drawArrow([0.95 0.5], 0.3, 1, 5, [0.8 0.4 0.4], 0.5);
        set([hl hr], 'visible', 'off');
    end

    function dialogResize(d, ~)
        [charsToPx, sizePx] = getCharsToPx(d);
        
        figEdgePx = round(data.settings.figEdge * charsToPx(1));
        factor = min((sizePx - figEdgePx) ./ fliplr(size(crtImage)));
        
        axesSize = round(fliplr(size(crtImage))*factor);
        
        set(ax, 'position', [0.5*(sizePx(1:2) - axesSize(1:2)) ...
            axesSize]);
    end

end