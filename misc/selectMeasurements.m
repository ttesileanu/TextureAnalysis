function measurements = selectMeasurements(measurements, mask)
% selectMeasurements Select a subset of a set of measurements.
%   sub = selectMeasurements(measurements, mask) uses the given mask to
%   select a subset of a measurements structure. The mask can be
%   index-based or logical, and is applied to each field of the
%   measurements structure.

fields = fieldnames(measurements);
for i = 1:length(fields)
    if isvector(measurements.(fields{i}))
        measurements.(fields{i}) = measurements.(fields{i})(mask);
    else
        measurements.(fields{i}) = measurements.(fields{i})(mask, :);
    end
end

end