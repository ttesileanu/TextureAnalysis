function color = mixcolor(color1, color2)
% Mix the two colors by adding them.
%   color = mixcolor(color1, color2) returns a color that is the sum of the
%   arguments, clipped to the range [0, 1].

color = max(min(color1 + color2, 1), 0);

end