function [mapping, details] = generateTextureMapping(fct, groups, directions, varargin)
% generateTextureMapping Generate mapping from one texture space to another.
%   mapping = generateTextureMapping(fct, groups, directions) generates
%   textures in the texture groups and directions identified by the cell
%   arrays `groups` and `directions`, then analyzes them using the function
%   `fct`, which should take an image patch and return a vector of texture
%   statistics. An interpolation is used to find the mapping from
%   coordinates along the directions of the input texture space to
%   coordinates in the target space (returned by `fct`). See the options to
%   control the number of samples, patch size, etc.
%
%   [mapping, details] = generateTextureMapping(...) also returns details
%   of the mapping, including the patches that were generated, and their
%   locations along the texture rays.
%
%   Options:
%    'nlocs'
%       Number of locations to generate along each direction.
%    'nsamples'
%       Number of texture samples to generate at each location.
%    'patchsize'
%       Size of the square texture patches that will be generated.
%    'order'
%       Interpolation order to use to generate the mapping. See
%       mapInterpolate.
%   
%   See also: mapInterpolate.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('nlocs', 8, @(n) isnumeric(n) && isscalar(n) && n >= 1);
parser.addParameter('patchsize', 64, @(n) isnumeric(n) && isscalar(n) && n >= 2);
parser.addParameter('nsamples', 16, @(n) isnumeric(n) && isscalar(n) && n >= 1);
parser.addParameter('order', 2, @(n) isnumeric(n) && isscalar(n) && n >= 1);

if nargin == 1 && strcmp(groups, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% setup the detailed output structure
n = length(groups);
details.locs = cell(n, 1);
details.patches = cell(n, 1);
details.stats = cell(n, 1);

% start a progress bar
progress = TextProgress('generating patches', 'prespace', 24, 'length', 20);
for i = 1:length(groups)
    % generate patches in this direction
    generator = PatchAxisGenerator(groups{i}, directions{i}, params.patchsize);
    generator.nLocations = params.nlocs;
    
    crt_patches = cell(1, params.nlocs);
    j = 1;
    while generator.next
        crt_patches{j} = generator.samples(params.nsamples);
        j = j + 1;
    end
    % store patches
    details.patches{i} = crt_patches;
    
    % analyze patches using target stats
    crt_evs = cell(1, params.nlocs);
    for j = 1:params.nlocs
        loc_ev = [];
        for p = 1:params.nsamples
            crt_patch = crt_patches{j}(:, :, p);
            crt_ev = fct(crt_patch);
            if isempty(loc_ev)
                loc_ev = zeros(params.nlocs, length(crt_ev));
            end
            loc_ev(p, :) = crt_ev; %#ok<AGROW>
        end
        crt_evs{j} = loc_ev;
    end
    % store stats
    details.stats{i} = crt_evs;
    % store locations
    details.locs{i} = generator.getLocations;
    
    progress.update(100*i/n);
end
progress.done('done');

% now generate mapping using interpolation
mapping = mapInterpolate(details.stats, details.locs, params.order);

end