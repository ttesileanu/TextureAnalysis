function [results, rawData] = loadTernaryPP(filename, varargin)
% loadTernaryPP Load ternary psychophysics results.
%   results = loadTernaryPP(filename) loads ternary psychophysics results
%   from the given file. By default, the function uses across-subjects
%   averages in the return values.
%
%   [results, rawData] = loadTernaryPP(...) also returns the raw contents
%   of the psychophysics file.
%
%   The results structure contains the following fields:
%    'groups': [n_meas, 1] cell array of strings
%       Cell array of coordinate groups for the thresholds. Multi-group
%       measurements are represented by separating the two group names by a
%       semicolon, while the direction in each group is identified by a
%       digit in parentheses. E.g., 'AB_1_1[0];AB_1_2[1]' is the mixed
%       group in which `p(A + B = 0 mod 3)` and `p(A + 2B = 1 mod 3)` are
%       varied.
%    'directions': [n_meas, 1] cell array of vectors
%       Cell array of normalized axis directions in which the thresholds
%       are measured. For single-group measurements, the axis is a
%       unit-norm three-component quasi-probability vector (i.e., summing
%       up to 1, but allowing negative components). For multi-group
%       measurements, this is a 3*ngroup-component vector in which each
%       group of 3 consecutive components adds up to 1. The normalization
%       is that from TERNARYDEC.
%    'thresholds': [n_meas, 1] vector
%       The threshold measurements, with 0 mapping to the fully-uncorrelated
%       center (i.e., all probability values equal to 1/3), and 1 mapping
%       to the normalized directions from the 'directions' array. This can
%       be NaN for missing measurements.
%    'thresholdIntervals': [n_meas, 2] matrix
%       Low (first column) and high (second column) 65% CI for the
%       thresholds. These can be NaN or infinite.
%    'nSubjects': [n_meas, 1] vector
%       Number of subjects that contribute for each measurement.
%    'subjects': [n_meas, 1] cell array
%       When the 'subjects' option is given, each direction can appear
%       several times, one for each subject for which measurements exist.
%       In that case, the 'subjects' field of the results structure
%       indicates which subject the respective measurement refers to.
%    'multi': [n_meas, 1] boolean vector
%       True for multi-group measurements.
%
%   The function also accepts a number of options:
%    'multi'
%       Set to 'true' to include multi-group planes in the output, 'false'
%       to focus on single-group planes.
%    'subjects'
%       Set to a string or cell array of strings to indicate which subjects
%       should be included in the output. The default, 'avg', returns only
%       average values. You can also use '*' to return all subjects.
%    'keepNaN'
%       Set to false to discard measurements for which the thresholds are
%       NaN.

% setup the input options
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

checkStr = @(s) isvector(s) && ischar(s);

parser.addParameter('multi', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('subjects', 'avg', @(s) checkStr(s) || (iscell(s) && ...
    all(cellfun(checkStr, s))));
parser.addParameter('keepNaN', true, @(b) isscalar(b) && islogical(b));

% show defaults
% XXX this can backfire if the user wants to process a file whose name is
% 'defaults'... hopefully that won't happen
if nargin == 1 && strcmp(filename, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% load the raw data
rawData = open(filename);

% expand options that can be both single strings and cell arrays
if ~iscell(params.subjects)
    if strcmp(params.subjects, '*')
        params.subjects = setdiff(fieldnames(rawData.ds_merged), {'avg'});
    else
        params.subjects = {params.subjects};
    end
end

% everything we will use is in ds_merged
data = rawData.ds_merged;

% start building the output structure
results.groups = {};
results.directions = {};
results.thresholds = [];
results.thresholdIntervals = [];
results.nSubjects = [];
results.subjects = {};
results.multi = [];

% go through the subjects one by one
for k = 1:length(params.subjects)
    crtSubject = params.subjects{k};
    
    % find all available texture groups for this subject
    subjData = data.(crtSubject);    
    subjEdirs = fieldnames(subjData.edirs);
    
    % mask out some unwanted/duplicate groups
    edirMask = true(size(subjEdirs));
    for i = 1:length(subjEdirs)
        edirData = subjData.edirs.(subjEdirs{i});
        cGroupNames = edirData.cgroup_names;
        
        % handle multi-group and single-group data differently
        if strcmp(cGroupNames{1}, cGroupNames{2})
            cGroupNames = cGroupNames{1};
            crtMulti = false;
        else
            crtMulti = true;
        end
        
        % skip multi-group data if necessary
        if ~params.multi && crtMulti
            edirMask(i) = false;
            continue;
        end
        
        % clean-up: if we have a choice, keep only fields that end in 'full'
        if ismember([subjEdirs{i} 'full'], subjEdirs)
            edirMask(i) = false;
            continue;
        end
        
        % clean-up: for A1 plane, keep only 'A1Gmerge'
        if ismember(subjEdirs{i}, {'A1G', 'A1Ginter'})
            % keep only merge
            edirMask(i) = false;
            continue;
        end
        
        % the suffix is different for the 'avg' subject
        if strcmp(crtSubject, 'avg')
            % XXX should ideally make configurable which of these we use
            crtSubjSuffix = '_all';
            crtNSubj = edirData.nsubjs_avail;
        else
            crtSubjSuffix = '';
            crtNSubj = 1;
        end
        
        % convert the 2-component data to 3-component axes for single
        % groups, or groups of 3-component axes for multi-group
        crtVecs = edirData.(['thresh_vecs' crtSubjSuffix]);
        crtVecsLo = edirData.(['thresh_vecs_eblo' crtSubjSuffix]);
        crtVecsHi = edirData.(['thresh_vecs_ebhi' crtSubjSuffix]);
        
        if ~crtMulti
            % single-group
            
            % convert to 3 components
            crtDirections = ternary2to3legacy(crtVecs);
            crtDirectionsLo = ternary2to3legacy(crtVecsLo);
            crtDirectionsHi = ternary2to3legacy(crtVecsHi);
            crtUVecs = ternary2to3legacy(edirData.uvecs);
        else
            % multi-group
            crtCGroupDirs = edirData.cgroup_dirs;
            
            crtDirections = ternarymix2to6(crtVecs, crtCGroupDirs);
            crtDirectionsLo = ternarymix2to6(crtVecsLo, crtCGroupDirs);
            crtDirectionsHi = ternarymix2to6(crtVecsHi, crtCGroupDirs);
            crtUVecs = ternarymix2to6(edirData.uvecs, crtCGroupDirs);            
        end
            
        % normalize and extract thresholds
        [crtThresholds, crtNormalized0] = ternarydec(crtDirections);
        [crtThreshLo, crtNormalizedLo] = ternarydec(crtDirectionsLo);
        [crtThreshHi, crtNormalizedHi] = ternarydec(crtDirectionsHi);
        
        % ensure that the normalized directions are the same
        tol = 1e-10;
        if any(abs(crtNormalized0(:) - crtNormalizedLo(:)) > tol)
            error([mfilename ':baddata'], 'Mismatch between normalized directions in threshold vs. eblo vectors.');
        end
        if any(abs(crtNormalized0(:) - crtNormalizedHi(:)) > tol)
            error([mfilename ':baddata'], 'Mismatch between normalized directions in threshold vs. ebhi vectors.');
        end
        
        % but prefer directions from uvecs, because those aren't
        % invalid in places where the thresholds are NaNs
        [~, crtNormalized] = ternarydec(crtUVecs);
        if any(abs(crtNormalized(:) - crtNormalized0(:)) > tol)
            error([mfilename ':baddata'], 'Mismatch between normalized directions in threshold vs. uvec vectors.');
        end
        
        if ~params.keepNaN
            % eliminate those measurements for which the threshold is NaN
            perDirMask = ~isnan(crtThresholds);
            if sum(perDirMask) == 0
                % completely skip subject, group pairs for which there
                % are no valid measurements
                continue;
            end
            crtThresholds = crtThresholds(perDirMask);
            crtThreshLo = crtThreshLo(perDirMask);
            crtThreshHi = crtThreshHi(perDirMask);
            crtNormalized = crtNormalized(perDirMask, :);
        end
        
        % turn directions into cell array
        crtNDirs = size(crtNormalized, 1);
        crtNormalizedCell = num2cell(crtNormalized, 2);
        
        if crtMulti
            % convert groups to single string
            crtGroupStr = [cGroupNames{1} '[' int2str(crtCGroupDirs(1)) '];' ...
                cGroupNames{2} '[' int2str(crtCGroupDirs(2)) ']'];
        else
            crtGroupStr = cGroupNames;
        end
        
        % add everything to the results
        results.groups = [results.groups ; repmat({crtGroupStr}, crtNDirs, 1)];
        results.directions = [results.directions ; crtNormalizedCell(:)];
        results.thresholds = [results.thresholds ; crtThresholds(:)];
        crtThreshIntervals = [crtThreshLo crtThreshHi];
        results.thresholdIntervals = [results.threshold_intervals ; crtThreshIntervals];
        results.subjects = [results.subjects ; repmat({crtSubject}, crtNDirs, 1)];
        results.multi = [results.multi ; repmat(crtMulti, crtNDirs, 1)];
        results.nSubjects = [results.n_subjects ; repmat(crtNSubj, crtNDirs, 1)];            
    end
end

% make sure the 'multi' array is logical (repmat makes it numeric)
results.multi = logical(results.multi);

end