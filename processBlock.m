function [P, ev, entropy, pattern_ix] = processBlock(I, nLevels)
% processBlock Calculate texture statistics from an image patch.
%   [P, ev, entropy, pattern_ix] = processBlock(I) calculates binary texture
%   statistics from an image patch by counting the occurrences of each
%   configuration of binary pixels for a square 2x2 glider. The 16 possible
%   configurations are then used to calculate 10 independent parameters
%   (see getStatistics).
%
%   processBlock(I, nLevels) calculates the statistics for nLevels gray
%   levels. XXX This isn't yet implemented.
%
%   Returns:
%    P:
%       Probabilities of gliders of each type. The encoding is as follows:
%           0               00             8               00
%                           00                             01
%
%           1               10             9               10
%                           00                             01 
%
%           2               01             10              01
%                           00                             01  
%
%           3               11             11              11
%                           00                             01  
%
%           4               00             12              00
%                           10                             11 
%
%           5               10             13              10
%                           10                             11  
%
%           6               01             14              01
%                           10                             11 
%
%           7               11             15              11
%                           10                             11
%    ev:
%       Probabilities projected onto the 10-dimensional space of
%       independent parameters (see getStatistics).
%    entropy:
%       Entropy of the discrete probability distribution identified by P.
%       (in bits)
%    pattern_ix:
%       Image in terms of the glider configurations. pattern_ix(i, j)
%       identifies the glider at I(i:i+1, j:j+1) using the code described
%       above.
%
%   See also: getStatistics.

if nargin < 2
    nLevels = 2;
end
if nLevels ~= 2
    error([mfilename ':notimp'], 'Only binary images are implemented for now.');
end

patterns = cat(3, ...
                   I(1:end-1,1:end-1), ...  % Upper left 
               2 * I(1:end-1,2:end), ...    % Upper right
               4 * I(2:end,1:end-1), ...    % Lower left
               8 * I(2:end,2:end));         % Lower right
               
pattern_ix = sum(patterns,3);

% The code:
%   0               00             8               00
%                   00                             01

%   1               10             9               10
%                   00                             01 

%   2               01             10              01
%                   00                             01  

%   3               11             11              11
%                   00                             01  

%   4               00             12              00
%                   10                             11 

%   5               10             13              10
%                   10                             11  

%   6               01             14              01
%                   10                             11 

%   7               11             15              11
%                   10                             11


% Restrict to non-overlapping patches (to avoid overcounting pixels)
%pattern_ix = pattern_ix(1:2:end,1:2:end);


P = histc(pattern_ix(:), 0:15);
P = P/sum(P);
%bar(0:15,P); pause(0.1);

%figure(2); bar(0:15,P); drawnow;
% Define the statistics:


stats = getStatistics;
ev = stats * P;
ix = find(P>0);
entropy = -sum(P(ix).*log2(P(ix)));


end