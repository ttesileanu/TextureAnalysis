function [stats,cons] = getStatistics(nLevels)
% getStatistics Generate mapping from all configurations of 2x2 gliders to
% independent components.
%   stats = getStatistics() generates a 10x16 matrix mapping all 16
%   possible configurations of a 2x2 binary glider to the 10 independent
%   components.
%
%   [stats, cons] = getStatistics() also returns a 6x16 matrix mapping the
%   16 glider probabilities to the 6 constraints. Projecting the
%   probabilities onto this 6-dimensional space should yield zero if the
%   probabilities were obtained from a valid image block.
%
%   The order for the 10 independent components is:
%    1. gamma
%    2. beta |
%    3. beta --
%    4. beta \
%    5. beta /
%    6. theta |-
%    7. theta _|
%    8. theta -|
%    9. theta |_
%   10. alpha
%   (see J. D. Victor and M. M. Conte,
%        "Local image statistics: maximum-entropy constructions and
%         perceptual salience,"
%        J. Opt. Soc. Am. A, vol. 29, no. 7, pp. 1313-1345, 2012)
%
%   The code for the 16 glider configuration is:
%       0               00             8               00
%                       00                             01
%
%       1               10             9               10
%                       00                             01 
%
%       2               01             10              01
%                       00                             01  
%
%       3               11             11              11
%                       00                             01  
%
%       4               00             12              00
%                       10                             11 
%
%       5               10             13              10
%                       10                             11  
%
%       6               01             14              01
%                       10                             11 
%
%       7               11             15              11
%                       10                             11
%
%   See also: processBlock.

if nargin < 1
    nLevels = 2;
end

if isfinite(nLevels) && nLevels ~= 2
    error([mfilename ':notimp'], 'Only binary and arbitrary grayscale images are implemented for now.');
end

stats = [...
        -1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1 ; ... % gamma
         1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1 ; ... % beta |
         1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1 ; ... % beta --
         1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1 ; ... %  beta \
         1, 1,-1,-1, -1,-1,1, 1, 1, 1,-1,-1,-1,-1, 1, 1 ; ... %  beta /
        -1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1 ; ... % theta |-
        -1,-1, 1, 1, 1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1 ; ... % theta _|
        -1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1 ; ... % theta -|
        -1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1 ; ... % theta |_
         1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1   ... % alpha
];

if nargout > 1
    cons = zeros(6,16);
    % Constraints
    
    % Equality of single-pixel marginals:
    % 0 + 2 + 4 + 6 + 8 + 10 + 12 + 14
    % - (0 + 1 + 4 + 5 + 8 + 9 + 12 + 13) = 0
    cons(1,:) = [0,-1,1,0,0,-1,1,0,0,-1,1,0,0,-1,1,0];
    
    
    %    0 +     2 +     4 + 6 + 8 +     10 +    12 + 14
    % - (0 + 1 + 2 + 3 +         8 + 9 + 10 + 11) = 0
    cons(2,:) = [0,-1,0,-1,1,0,1,0,0,-1,0,-1,1,0,1,0];
    
    %    0 +     2 +     4 +     6 +    8 +     10 +    12 + 14
    % -( 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7)
    cons(3,:) = [0,-1,0,-1,0,-1,0,-1,1,0,1,0,1,0,1,0];
    
    % Equality of two-pixel marginals:
    % Vertical:
    %  (0 +     2 +       8 + 10)
    %- (0 + 1 +     4 + 5) = 0
    cons(4,:) = [0,-1,1,0,-1,-1,0,0,1,0,1,0,0,0,0,0];
    
    % Horizontal:
    % 0 + 4 + 8 + 12
    % - (0 + 1 + 2 + 3)
    cons(5,:) = [0,-1,-1,-1,1,0,0,0,1,0,0,0,1,0,0,0];
    
    
    cons(6,:) = 1;                                                              % Normalization
end

end