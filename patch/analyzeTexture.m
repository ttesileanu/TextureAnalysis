function [ev, P, entropy, patternIdx] = analyzeTexture(I, nLevels)
% analyzeTexture Calculate texture statistics from an image patch.
%   [ev, P, entropy, patternIdx] = analyzeTexture(I, 2) calculates binary
%   texture statistics from an image patch by counting the occurrences of
%   each configuration of binary pixels for a square 2x2 glider. The 16
%   possible configurations, returned in the vector `P`, are then used to
%   calculate 10 independent parameters, returned in `ev`. The entropy of
%   the pixels in the patch is also returend, together with `patternIdx`,
%   which identifies the block at each location, i.e., `P` is a histogram
%   of `patternIdx` values.
%
%   The ordering of the texture coordinates in `ev` is the following:
%      1. gamma    ; 2. beta |   ; 3. beta --  ; 4. beta \   ;  5. beta /
%      6. theta |- ; 7. theta _| ; 8. theta -| ; 9. theta |_ ; 10. alpha
%   Each coordinate corresponds to a glider shape. For each shape, we count
%   how many times the glider overlaps an even number of pixels equal to 1,
%   nEven, and how many times it overlaps an odd number, nOdd, as it slides
%   over the entire patch. The value of the corresponding binary texture
%   coordinate is given by (nEven - nOdd) / (nEven + nOdd).
%
%   The mapping from 2x2 blocks to values in `patternIdx` (and thus indices
%   in `P`) is given by multiplying the blocks element by element with
%                       [1 2]
%                       [4 8]
%   and then summing the resulting values.
%
%   [ev, P] = analyzeTexture(I, nLevels) calculates texture statistics for
%   the patch `I` using the graylevel analysis with `nLevels` levels. The
%   probability vector of all `nLevels^4` configurations of 2x2 blocks of
%   pixels is returned in `P`, while `ev` contains a representation in
%   terms of independent coordinates. Specifically, the dimensions go
%   through the texture groups in the order from `mtc.coord_groups`, with
%   `nLevels - 1` values for each group representing all the probability
%   values except for the last 1 (which is redundant, since they add up to
%   1). (The `mtc` structure here can be obtained by running
%   `processBlock('mtc', nLevels).) Note that this representation is
%   different from that used in the binary case (which is also used for the
%   continuous textures described below)!
%
%   analyzeTexture(I, inf) calculates continuous texture statistics. These
%   are combinations of 1-, 2-, 3-, and 4-point correlation functions
%   designed to exactly recover the binary texture results when the input
%   `I` is binary (entries 0 or 1). One way to think about the texture
%   coordinates in this case is to imagine stochastically binarizing the
%   grayscale image by using the pixel value `I(i, j)` as the probability
%   of setting that pixel to white (=1). If we average the binary texture
%   coordinates obtained for each such binarized patch over the ensemble of
%   all binarizations, we get the continuous statistics.
%
%   Note that if `nLevels` is different from 2 or inf, functions from
%   Jonathan Victor's framework are used and required to be on the path
%   (specifically, these are gtc_define, mtc_define, glider_mapubi, and
%   mtc_probs2cgs).
%
%   References:
%     Hermundstad et al. (2014). Variance predicts salience in central sensory processing. ELife, 10.7554, e03722.
%     Victor et al. (2012). Local image statistics: maximum-entropy constructions and perceptual salience. Journal of the Optical Society of America. A, Optics, Image Science, and Vision, 29(7), 1313?45.

% set up the persistent structures for Jonathan's code
% these are used in conjunction with Jonahtan's code when nLevels isn't inf or 2
persistent gtc;
persistent mtcs;

% figure out whether we need to generate gtc or mtc structures
need_gtc = false;
need_mtc = false;

% are we explicitly asking for gtc and mtc structures?
if strcmp(I, 'gtc')
    need_gtc = true;
elseif strcmp(I, 'mtc')
    % we need to generate the gtc in order to get mtc
    need_gtc = true;
    need_mtc = true;
elseif isfinite(nLevels) && nLevels > 2
    % we need gtc and mtc if we're working with images with finite number
    % of gray levels that's larger than 2
    need_gtc = true;
    need_mtc = true;
end

% generate gtc if needed
checks = [0 0 ; 0 1 ; 1 0 ; 1 1]; % used also in glider_mapubi below
if need_gtc && isempty(gtc)
    gtc = gtc_define(checks);
end

% generate mtc if needed
if need_mtc
    if length(mtcs) < nLevels
        mtcs{nLevels} = [];
    end
    if isempty(mtcs{nLevels})
        mtcs{nLevels} = mtc_define(nLevels, gtc);
    end
end

% if this is a request for gtc or mtc, answer it
if ischar(I) && isvector(I)
    switch I
        case 'gtc'
            ev = gtc;
            return;
        case 'mtc'
            ev = mtcs{nLevels};
            return
        otherwise
            error([mfilename ':badreq'], 'Unknown request.');
    end
end

% if this is a request to process an image patch, do it
if isfinite(nLevels) && nLevels > 2
    % calculate stats for current patch
    P = glider_mapubi(round(I*(nLevels-1)), checks, nLevels);
    P = flatten(P/sum(P));
    
    % calculate independent stats
    evmat = mtc_probs2cgs(P, mtcs{nLevels});
    % remove one of the columns, since it is always equal to 1 minus the others
    ev = flatten(evmat(:, 1:end-1)');
    
    return;
end

% use the simpler code for 2 gray levels, and the continuous stats code for
% nLevels == inf
if nLevels == 2
    patterns = cat(3, ...
        I(1:end-1, 1:end-1), ...        % upper left
        2 * I(1:end-1, 2:end), ...      % upper right
        4 * I(2:end, 1:end-1), ...      % lower left
        8 * I(2:end, 2:end));           % lower right
    
    patternIdx = sum(patterns, 3);
    
    P = histcounts(patternIdx(:), 0:16);
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

stats = [...
        -1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1 ; ... % gamma
         1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1 ; ... % beta |
         1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1 ; ... % beta --
         1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1 ; ... % beta  \
         1, 1,-1,-1, -1,-1,1, 1, 1, 1,-1,-1,-1,-1, 1, 1 ; ... % beta  /
        -1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1 ; ... % theta |-
        -1,-1, 1, 1, 1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1 ; ... % theta _|
        -1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1 ; ... % theta -|
        -1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1 ; ... % theta |_
         1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1   ... % alpha
];

ev = P * stats';

if nLevels == 2 && nargin > 2
    Pix = P(P>0);
    entropy = -sum(Pix.*log2(Pix));
end

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