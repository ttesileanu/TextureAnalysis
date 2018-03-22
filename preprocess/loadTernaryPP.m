function results = loadTernaryPP(filename, varargin)
% loadTernaryPP Load ternary psychophysics results.
%   results = loadTernaryPP(filename) loads ternary psychophysics results
%   from the given file. By default, the function uses across-subjects
%   averages in the return values.
%
%   The results structure contains the following fields:
%    'groups': [n_meas, 1] cell array
%       Cell array of coordinate groups for the thresholds. For multi-group
%       measurements, each entry in this cell array will itself be a cell
%       array, listing all the groups involved.
%    'directions': [n_meas, 1] cell array
%       Cell array of normalized axis directions in which the thresholds
%       are measured. For single-group measurements, the axis is a
%       unit-norm three-component pseudo-probability vector (i.e., summing
%       up to 1, but allowing negative components). For multi-group
%       measurements, this is a cell array with one 3-component vector for
%       each group.
%    'thresholds': [n_meas, 1] vector
%       The threshold measurements, with 0 mapping to the
%       fully-uncorrelated center (i.e., all probability values equal to
%       1/3), and 1 mapping to the normalized directions from the
%       'directions' array. This can be NaN.
%    'threshold_intervals': [n_meas, 2] matrix
%       Low (first column) and high (second column) 65% CI for the
%       thresholds. These can be NaN or infinite.
%    'n_subjects': [n_meas, 1] vector
%       Number of subjects that contribute for each measurement.
%    'subjects': [n_meas, 1] cell array
%       When the 'subjects' option is given, each direction can appear
%       several times, one for each subject for which measurements exist.
%       In that case, the 'subjects' field of the results structure
%       indicates which subject the respective measurement refers to.
%
%   The function also accepts a number of options:
%    'multi'
%       Set to 'true' to include multi-group planes in the output, 'false'
%       to focus on single-group planes.
%    'subjects'
%       Set to a string or cell array of strings to indicate which subjects
%       should be included in the output. The default, 'all', returns only
%       average values. You can also use '*' to return all subjects.
%    'exclude'
%       Cell array of subjects to exclude.
%    'keepnan'
%       Set to false to discard measurements for which the thresholds are
%       NaN.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

checkStr = @(s) isvector(s) && ischar(s);

parser.addParameter('multi', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('subjects', 'all', @(s) checkStr(s) || (iscell(s) && ...
    all(cellfun(checkStr, s))));
parser.addParameter('exclude', {}, @(s) checkStr(s) || (iscell(s) && ...
    all(cellfun(checkStr, s))));
parser.addParameter('keepnan', true, @(b) isscalar(b) && islogical(b));

% XXX this can backfire if the user wants to process a file whose name is
% 'defaults'... hopefully that seems unlikely
if nargin == 1 && strcmp(filename, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% expand options that can be both single strings and cell arrays
if ~iscell(params.subjects)
    if isempty(params.subjects)
        params.subjects = {};
    else
        params.subjects = {params.subjects};
    end
end

if ~iscell(params.exclude)
    if isempty(params.exclude)
        params.exclude = {};
    else
        params.exclude = {params.exclude};
    end
end

% load the raw data
raw_data = open(filename);

% using the 'ds_merged.avg' data for everything, since it also contains
% per-subject information
data = raw_data.ds_merged.avg;

% find all the texture groups
all_edirs = fieldnames(data.edirs);

% mask out unwanted texture groups
mask = true(size(all_edirs));
for i = 1:length(all_edirs)
    crt_edir_data = data.edirs.(all_edirs{i});
    crt_cgroup_names = crt_edir_data.cgroup_names;
    % skip multi-group data if necessary
    if ~params.multi && ~strcmp(crt_cgroup_names{1}, crt_cgroup_names{2})
        mask(i) = false;
        continue;
    end
    
    % clean-up: if we have a choice, keep only fields that end in 'full'
    if ismember([all_edirs{i} 'full'], all_edirs)
        mask(i) = false;
        continue;
    end

    % clean-up: for A1 plane, keep only 'A1Gmerge'
    if ismember(all_edirs{i}, {'A1G', 'A1Ginter'})
        % keep only merge
        mask(i) = false;
        continue;
    end
end
sub_edirs = all_edirs(mask);

% start building the output structure
results.groups = {};
results.directions = {};
results.thresholds = [];
results.threshold_intervals = [];
results.n_subjects = [];
results.subjects = {};

% process the raw data
for i = 1:length(sub_edirs)
    crt_edir_data = data.edirs.(sub_edirs{i});
    crt_groups = crt_edir_data.cgroup_names;
    
    % handle multi-group and single-group data differently, for consistency
    % with older code
    if strcmp(crt_groups{1}, crt_groups{2})
        crt_groups = crt_groups{1};
    end
    
    % figure out which subjects to include
    if length(params.subjects) == 1 && ismember(params.subjects{1}, {'all', '*'})
        if ~strcmp(params.subjects{1}, '*')
            crt_subjects = params.subjects;
        else
            crt_subjects = crt_edir_data.subjs_avail_ids;
        end
    else
        crt_subjects = intersect(crt_edir_data.subjs_avail_ids, params.subjects);
    end
    
    % figure out which subjects to exclude
    crt_subjects = setdiff(crt_subjects, params.exclude);
    
    % go through the subjects and process the data
    for j = 1:length(crt_subjects)
        crt_subject = crt_subjects{j};
        
        % convert the 2-component data to 3-component axes for single
        % groups, or groups of 3-component axes for multi-group
        crt_vecs = crt_edir_data.(['thresh_vecs_' crt_subject]);
        crt_vecs_lo = crt_edir_data.(['thresh_vecs_eblo_' crt_subject]);
        crt_vecs_hi = crt_edir_data.(['thresh_vecs_ebhi_' crt_subject]);
        
        if ischar(crt_groups)
            % single-group
            
            % convert to 3 components
            crt_directions = ternary2to3legacy(crt_vecs);
            crt_directions_lo = ternary2to3legacy(crt_vecs_lo);
            crt_directions_hi = ternary2to3legacy(crt_vecs_hi);
            
            % normalize and extract thresholds
            [crt_thresh, crt_normalized0] = ternarydec(crt_directions);
            [crt_thresh_lo, crt_normalized_lo] = ternarydec(crt_directions_lo);
            [crt_thresh_hi, crt_normalized_hi] = ternarydec(crt_directions_hi);
            
            % ensure that the normalized directions are the same
            tol = 1e-10;
            if any(abs(crt_normalized0(:) - crt_normalized_lo(:)) > tol)
                error([mfilename ':baddata'], 'Mismatch between normalized directions in threshold vs. eblo vectors.');
            end
            if any(abs(crt_normalized0(:) - crt_normalized_hi(:)) > tol)
                error([mfilename ':baddata'], 'Mismatch between normalized directions in threshold vs. ebhi vectors.');
            end
            
            % but prefer directions from uvecs, because those aren't
            % invalid in places where the thresholds are NaNs
            [~, crt_normalized] = ternarydec(ternary2to3legacy(crt_edir_data.uvecs));
            if any(abs(crt_normalized(:) - crt_normalized0(:)) > tol)
                error([mfilename ':baddata'], 'Mismatch between normalized directions in threshold vs. uvec vectors.');
            end
            
            if ~params.keepnan
                % eliminate those measurements for which the threshold is NaN
                mask = ~isnan(crt_thresh);
                if sum(mask) == 0
                    % completely skip subject, group pairs for which there
                    % are no valid measurements
                    continue;
                end
                crt_thresh = crt_thresh(mask);
                crt_thresh_lo = crt_thresh_lo(mask);
                crt_thresh_hi = crt_thresh_hi(mask);
                crt_normalized = crt_normalized(mask, :);
            end
            
            % turn directions into cell array
            crt_ndirs = size(crt_normalized, 1);
            crt_normalized_cell = mat2cell(crt_normalized, ones(crt_ndirs, 1), 3);
            
            % add everything to the results
            results.groups = [results.groups ; repmat({crt_groups}, crt_ndirs, 1)];
            results.directions = [results.directions ; crt_normalized_cell(:)];
            results.thresholds = [results.thresholds ; crt_thresh(:)];
            crt_thresh_intervals = [crt_thresh_lo crt_thresh_hi];
            results.threshold_intervals = [results.threshold_intervals ; crt_thresh_intervals];
            results.subjects = [results.subjects ; repmat({crt_subject}, crt_ndirs, 1)];
            if strcmp(crt_subject, 'all')
                crt_n_subj = crt_edir_data.nsubjs_avail;
            else
                crt_n_subj = 1;
            end
            results.n_subjects = [results.n_subjects ; repmat(crt_n_subj, crt_ndirs, 1)];
        else
            % multi-group
            % XXX don't know exactly which directions are used in the
            % multi-group cases
            error([mfilename ':notimp'], 'Multi-group processing not yet implemented.');
        end
    end
end

end