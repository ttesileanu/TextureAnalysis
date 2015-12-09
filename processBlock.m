function [P,ev,entropy,pattern_ix] = processBlock(I)

% Collect histogram of patterns for square glider

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

  
  