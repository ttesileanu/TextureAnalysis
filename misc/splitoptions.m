function varargout = splitoptions(options, varargin)
% SPLITOPTIONS Split a cell array of options in the format string, value
% into groups based on the option names.
%   [opts1, ..., optsN, others] = SPLITOPTIONS(options, names1, ..., namesN)
%   splits the cell array of options into several groups. The options in the
%   cell array must be pairs in the format string, value. All of the options
%   whose strings match one of the names in namesX will be placed in optsX.
%   Overlaps are allowed between the names in the different groups. Options
%   whose names appear in none of the lists are placed in the 'others'
%   group.

if ~isempty(options) && (~iscell(options) || ~isvector(options) || mod(length(options), 2) ~= 0 || ...
        ~all(arrayfun(@(i) ischar(options{i}) && isvector(options{i}), 1:2:length(options))))
    error([mfilename ':badoptions'], 'The options should come in pairs of the form string, value.');
end

varargout = cell(length(varargin) + 1, 1);
allmasks = false(size(options));
allmasks = allmasks(:);
for i = 1:length(varargin)
    if ~all(cellfun(@(s) ischar(s) && isvector(s), varargin{i}))
        error([mfilename ':badnames'], 'The names should be cell arrays of strings.');
    end
    idxs = find(ismember(options(1:2:end), varargin{i}));
    mask = false(size(allmasks));
    mask(2*idxs-1) = true;
    % include the values
    mask(2*idxs) = true;
    varargout{i} = options(mask);
    allmasks = (allmasks | mask);
end
varargout{length(varargin) + 1} = options(~allmasks);

end