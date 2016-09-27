function c = structToCell(s, sel)
% structToCell Create cell array of key-value pairs from structure.
% c = structToCell(s) creates a cell array of the form {key1, value1, key2,
% value2, ...} representing the fields of the structure s.
% 
% structToCell(s, sel) focuses on only the fields whose keys match the
% selection.
%
% Fields with empty values are ignored.

fields = fieldnames(s);
if nargin > 1
    fields = intersect(fields, sel);
end

mask = cellfun(@(field) ~isempty(s.(field)), fields);
fields = fields(mask);

c = cell(1, 2*length(fields));
c(1:2:end) = fields;
c(2:2:end) = cellfun(@(field) s.(field), fields, 'uniform', false);

end