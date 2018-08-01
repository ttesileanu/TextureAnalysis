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
        if isempty(measurements.(field))
            if isfield(varargin{i}, field)
                measurements.(field) = varargin{i}.(field);
            end
        else
            if isfield(varargin{i}, field)
                measurements.(field) = [measurements.(field) ; varargin{i}.(field)];
            else
                switch field
                    case {'thresholdIntervals', 'nSubjects'}
                        measurements.(field) = [measurements.(field) ; ...
                            nan(length(varargin{i}.groups), size(measurements.(field), 2))];
                    case 'subjects'
                        measurements.(field) = [measurements.(field) ; ...
                            repmat({''}, length(varargin{i}.groups), 1)];
                    case 'multi'
                        measurements.(field) = [measurements.(field) ; ...
                            cellfun(@(g) 1 + sum(g == ';'), varargin{i}.groups(:))];
                    otherwise
                        error([mfilename ':badfield'], 'Don''t know how to fill in default values for %s.', field);
                end
            end
        end
    end
end

end