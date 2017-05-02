function plotChange(x1, y1, x2, y2, varargin)
% plotChange Show changes from one scatter plot to another.
%   plotChange(x1, y1, x2, y2) shows the changes from the dataset (x1, y1)
%   to (x2, y2). x1, y1, x2, and y2 must be vectors of the same length.
%   (x1, y1) is in red, going towards blue for (x2, y2).

n = length(x1);
if ~isvector(x1) || ~isvector(y1) || ~isvector(x2) || ~isvector(y2) || ...
        length(x2) ~= n || length(y1) ~= n || length(y2) ~= n
    error([mfilename ':badinp'], 'All input data must be vectors of the same size.');
end

% XXX make configurable
col1 = [0.1 0.1 1.0];
col2 = [0.8 0.2 0.2];
alpha = 0.4;

cmap = zeros(64, 3);
for k = 1:64
    a = (k - 1)/63;
    cmap(k, :) = (1-a)*col1 + a*col2;
end

hold on;
colormap(cmap);
% nsegs = 3;
for i = 1:n
%     v1 = [x1(i) y1(i)];
%     v2 = [x2(i) y2(i)];
    
   crtX = [x1(i) x2(i)];
   crtY = [y1(i) y2(i)];
   crtZ = zeros(1, 2);
   crtCol = [0 1];
    
   surface([crtX; crtX], [crtY; crtY], [crtZ; crtZ], [crtCol; crtCol], ...
       'facecolor', 'none', 'edgecolor', 'interp', 'edgealpha', alpha);

%     for k = 1:nsegs
%         a = (k - 1) / nsegs;
%         b = k / nsegs;
%         w1 = (1 - a)*v1 + a*v2;
%         w2 = (1 - b)*v1 + b*v2;
%         col = cmap(1 + floor((size(cmap, 1) - 1)*(k - 1) / (nsegs - 1)), :);
%         
%         plot([w1(1) w2(1)], [w1(2) w2(2)], 'color', col);
%     end
end

h1 = plot(x1, y1, '.', 'color', col1);
h2 = plot(x2, y2, '.', 'color', col2);

drawnow;

% XXX undocumented Matlab, because they're too stupid to expose the
% transparency properties for LINE objects returned by plot...
h1Markers = h1.MarkerHandle;
h2Markers = h2.MarkerHandle;

h1Markers.EdgeColorData = uint8([h1Markers.EdgeColorData(1:3) ; round(alpha*255)]);
h2Markers.EdgeColorData = uint8([h2Markers.EdgeColorData(1:3) ; round(alpha*255)]);

end