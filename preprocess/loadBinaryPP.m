function [results, raw_data] = loadBinaryPP(filename, varargin)
% loadBinaryPP Load binary psychophysics results.
%   results = loadBinaryPP(filename) loads binary psychophysics results
%   from the given file. By default, the function uses across-subjects
%   averages in the return values.
%
%   [results, raw_data] = loadTernaryPP(...) also returns the raw contents
%   of the psychophysics file.
%
%   The results structure contains the following fields:
%    'groups': [n_meas, 1] cell array
%       Cell array of coordinate groups for the thresholds. For multi-group
%       measurements, each entry in this cell array will itself be a cell
%       array, listing all the groups involved.
%    'directions': [n_meas, 1] cell array
%       Cell array of normalized axis directions in which the thresholds
%       are measured. For single-group measurements, the axis is a
%       unit-norm two-component quasi-probability vector (i.e., summing
%       up to 1, but allowing negative components). For multi-group
%       measurements, this is a 2*ngroup-component vector in which each
%       group of 2 consecutive components adds up to 1. The normalization
%       is that from TERNARYDEC.
%    'thresholds': [n_meas, 1] vector
%       The threshold measurements, with 0 mapping to the
%       fully-uncorrelated center (i.e., both probability values equal to
%       1/2), and 1 mapping to the normalized directions from the
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
%    'multi': [n_meas, 1] boolean vector
%       True for multi-group measurements.
%
%   The function also accepts a number of options:
%    'multi'
%       Set to 'true' to include multi-group planes in the output, 'false'
%       to focus on single-group measurements.
%    'subjects'
%       Set to a string or cell array of strings to indicate which subjects
%       should be included in the output. The default, 'avg', returns only
%       average values. You can also use '*' to return all subjects.
%    'exclude'
%       Cell array of subjects to exclude.
%    'keepnan'
%       Set to false to discard measurements for which the thresholds are
%       NaN.
%    'ninterp'
%       Number of points where to sample the mixed-group psychophysics
%       ellipses within each quadrant. That is, if this is set to 1, each
%       plane will have 4 points. Sampling is done at equal angles in the
%       plane.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

checkStr = @(s) isvector(s) && ischar(s);

parser.addParameter('multi', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('subjects', 'avg', @(s) checkStr(s) || (iscell(s) && ...
    all(cellfun(checkStr, s))));
parser.addParameter('exclude', {}, @(s) checkStr(s) || (iscell(s) && ...
    all(cellfun(checkStr, s))));
parser.addParameter('keepnan', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('ninterp', 1, @(n) isscalar(n) && isnumeric(n) && n >= 1);

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

% load the raw data
raw_data = open(filename);
pp_data = raw_data.dataPP;

% expand options that can be both single strings and cell arrays
if ~iscell(params.subjects)
    if isempty(params.subjects)
        params.subjects = {};
    else
        if strcmp(params.subjects, '*')
            params.subjects = {pp_data.indS.subjID};
        else
            params.subjects = {params.subjects};
        end
    end
end

if ~iscell(params.exclude)
    if isempty(params.exclude)
        params.exclude = {};
    else
        params.exclude = {params.exclude};
    end
end

% start building the output structure
results.groups = {};
results.directions = {};
results.thresholds = [];
% results.threshold_intervals = [];
% results.n_subjects = [];
results.subjects = {};
results.multi = [];

% create the structure of directions we will be looking at

% single groups
% these are ordered in the same way as the coordinates in the raw data
edirs1 = {'A_1', 'AC_1_1', 'AB_1_1', 'AD_1_1', 'BC_1_1', ...
    'ABC_1_1_1', 'BCD_1_1_1', 'ABD_1_1_1', 'ACD_1_1_1', 'ABCD_1_1_1_1'};

% decide which angles in the mixed planes we use for sampling
rem_last = @(v) v(1:end-1);
mixed_angles = flatten(bsxfun(@plus, [0 pi/2 pi 3*pi/2], ...
    rem_last(linspace(0, pi/2, params.ninterp+1))'));
mixed_trig = [cos(mixed_angles') ; sin(mixed_angles')];

% process the raw data for all subjects
for k = 1:length(params.subjects)
    % collect covariance matrix for the appropriate subject
    crt_subject = params.subjects{k};
    if ismember(crt_subject, params.exclude)
        continue;
    end
    if strcmp(crt_subject, 'avg')
        % the factor 1000 was used in the encoding
        crt_edir_data = 1000*pp_data.allS.subjAvg.covM;
    else
        % the factor 1000 was used in the encoding
        subj_mask = strcmp({pp_data.indS.subjID}, crt_subject);
        crt_edir_data = 1000*pp_data.indS(subj_mask).covM;
    end
    
    % iterate through all group pairs
    for i = 1:length(edirs1)
        group1 = edirs1{i};
        for j = 1:i-1
            group2 = edirs1{j};
            
            % focus on the appropriate 2x2 component of the covariance matrix
            crt_cov2 = crt_edir_data([i j], [i j]);
            
            % generated interpolated locations
            m = sqrt(1 ./ sum(mixed_trig.* (crt_cov2*mixed_trig), 1));
            interpolated0 = bsxfun(@times, m, mixed_trig);
            
            % p_x = [(1-x)/2, (1+x)/2]
            % p_y = [(1-y)/2, (1+y)/2]
            % norm(p_x)^2 + norm(p_y)^2 = 1/4*(1 - 2x + x^2 + 1 + 2x + x^2 + 
            %   1 - 2y + y^2 + 1 + 2y + y^2) = 1 + (x^2 + y^2)/2
            
            % map locations to quasi-probability space
            interpolated = binary1to2(interpolated0');
            
            % convert to magnitude and direction
            [mags, uvecs] = binarydec(interpolated);
            
            % convert directions to cell array
            crt_directions = num2cell(uvecs, 2);
            
            % add everything to the results
            crt_ndirs = length(mags);
            crt_group_str = [group1 ';' group2];
            results.groups = [results.groups ; repmat({crt_group_str}, crt_ndirs, 1)];
            results.directions = [results.directions ; crt_directions(:)];
            results.thresholds = [results.thresholds ; mags(:)];
            results.subjects = [results.subjects ; repmat({crt_subject}, crt_ndirs, 1)];
            results.multi = [results.multi ; true(crt_ndirs, 1)];
        end
        
        % add single-group results
        % [v v ; v v]*[1 ; 0] = [v ; v] .* [1 ; 0] = [v ; 0]
        % sum([v ; 0], 1) = v
        % m = 1 ./ sqrt(v)
        % interp = [m ; 0] = [1./sqrt(v) ; 0]
        mags = 1./sqrt(crt_edir_data(i, i));
        crt_directions = [1 0];
        crt_group_str = group1;
        results.groups = [results.groups ; crt_group_str];
        results.directions = [results.directions ; crt_directions];
        results.thresholds = [results.thresholds ; mags(:)];
        results.subjects = [results.subjects ; {crt_subject}];
        results.multi = [results.multi ; false];
    end
end

end