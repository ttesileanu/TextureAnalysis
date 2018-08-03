function measurements = catMeasurements(varargin)
% catMeasurements Concatenate measurements structures.
%   measurements = catMeasurements(measurements1, measurements2, ...)
%   concatenates the given measurements structures. Fields that don't exist
%   in some of the structures are filled with default values: NaN for
%   `thresholdIntervals` and `nSubjects`, empty strings for `subjects`, and
%   a value inferred from the group name for `multi`.

% find out which fields we need to have in the output
allFields = {};
for i = 1:length(varargin)
    allFields = union(allFields, fieldnames(varargin{i}), 'stable');
end

% create the ouptut structure
measurements = struct;
for k = 1:length(allFields)
    field = allFields{k};
    measurements.(field) = [];
    for i = 1:nargin
        nCrt = length(varargin{i}.groups);
        if isfield(varargin{i}, field)
            value = varargin{i}.(field);
        else
            switch field
                case 'thresholdIntervals'
                    value = nan(nCrt, 2);
                case 'nSubjects'
                    value = nan(nCrt, 1);
                case 'subjects'
                    value = repmat({''}, nCrt, 1);
                case 'multi'
                    value = (cellfun(@(g) sum(g == ';'), varargin{i}.groups(:)) > 0);
                otherwise
                    error([mfilename ':badfield'], 'Don''t know how to fill in default values for %s.', field);
            end
        end
        
        if isempty(measurements.(field))
            measurements.(field) = value;
        else
            measurements.(field) = [measurements.(field) ; value];
        end
    end
end

end
