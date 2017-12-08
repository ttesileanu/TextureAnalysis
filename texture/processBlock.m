function [P, ev, entropy, pattern_ix] = processBlock(I, nLevels)
% processBlock Calculate texture statistics from an image patch.
%   [P, ev, entropy, pattern_ix] = processBlock(I) calculates binary texture
%   statistics from an image patch by counting the occurrences of each
%   configuration of binary pixels for a square 2x2 glider. The 16 possible
%   configurations are then used to calculate 10 independent parameters
%   (see getStatistics).
%
%   processBlock(I, nLevels) calculates the statistics for nLevels gray
%   levels. XXX This isn't yet implemented (unless nLevels is 2 or inf).
%
%   processBlock(I, inf) calculates the statistics for full grayscale
%   images. These are combinations of 1-, 2-, 3-, and 4-point correlation
%   functions designed to exactly recover the binary texture results when
%   the input I is binary (entries 0 or 1).
%
%   Returns:
%    P:
%       Probabilities of gliders of each type. The encoding for binary
%       statistics is as follows:
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
%       above. This is only defined when nLevels == 2.
%
%   See also: getStatistics.

if nargin < 2
    nLevels = 2;
end
if isfinite(nLevels) && nLevels ~= 2
    error([mfilename ':notimp'], 'Only binary and arbitrary grayscale images are implemented for now.');
end

if nLevels == 2
    patterns = cat(3, ...
        I(1:end-1,1:end-1), ...  % Upper left
        2 * I(1:end-1,2:end), ...    % Upper right
        4 * I(2:end,1:end-1), ...    % Lower left
        8 * I(2:end,2:end));         % Lower right
    
    pattern_ix = sum(patterns,3);
    
    P = histc(pattern_ix(:), 0:15);
    P = P/sum(P);
elseif ~isfinite(nLevels)
    % we will need shifted&flattened versions of the image patch
    I1 = flatten(I(1:end-1, 1:end-1));
    I2 = flatten(I(1:end-1, 2:end));
    I3 = flatten(I(2:end, 1:end-1));
    I4 = flatten(I(2:end, 2:end));
    
    % next, we will need means of all the combinations of products
    I1m = mean(I1);
    I2m = mean(I2);
    I3m = mean(I3);
    I4m = mean(I4);
    
    I12m = mean(I1 .* I2);
    I13m = mean(I1 .* I3);
    I14m = mean(I1 .* I4);
    I23m = mean(I2 .* I3);
    I24m = mean(I2 .* I4);
    I34m = mean(I3 .* I4);
    
    I123m = mean(I1 .* I2 .* I3);
    I124m = mean(I1 .* I2 .* I4);
    I134m = mean(I1 .* I3 .* I4);
    I234m = mean(I2 .* I3 .* I4);
    
    I1234m = mean(I1 .* I2 .* I3 .* I4);
    
    P = zeros(16, 1);
    P( 1) = 1 - I1m - I2m - I3m - I4m + I12m + I13m + I14m + I23m + I24m + I34m ...
                - I123m - I124m - I134m - I234m + I1234m;
    P( 2) = I1m - I12m - I13m - I14m + I123m + I124m + I134m - I1234m;
    P( 3) = I2m - I12m - I23m - I24m + I123m + I124m + I234m - I1234m;
    P( 4) = I12m - I123m - I124m + I1234m;
    P( 5) = I3m - I13m - I23m - I34m + I123m + I134m + I234m - I1234m;
    P( 6) = I13m - I123m - I134m + I1234m;
    P( 7) = I23m - I123m - I234m + I1234m;
    P( 8) = I123m - I1234m;
    P( 9) = I4m - I14m - I24m - I34m + I124m + I134m + I234m - I1234m;
    P(10) = I14m - I124m - I134m + I1234m;
    P(11) = I24m - I124m - I234m + I1234m;
    P(12) = I124m - I1234m;
    P(13) = I34m - I134m - I234m + I1234m;
    P(14) = I134m - I1234m;
    P(15) = I234m - I1234m;
    P(16) = I1234m;
end

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

%bar(0:15,P); pause(0.1);

%figure(2); bar(0:15,P); drawnow;
% Define the statistics:

if nargin > 1
    stats = getStatistics(nLevels);
    ev = stats * P;
    if nLevels == 2 && nargin > 2
        Pix = P(P>0);
        entropy = -sum(Pix.*log2(Pix));
    end
end

end