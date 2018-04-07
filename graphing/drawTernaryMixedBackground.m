function drawTernaryMixedBackground(groupdir1, groupdir2)
% drawTernaryMixedBackground Draw guiding axes through origin and guiding
% circles for pairs of ternary texture groups.
%
%   drawTernaryMixedBackground(groupdir1, groupdir2) draws guiding axes
%   going through the origin and guiding circles for a plot representing a
%   mix of two ternary texture groups. The arguments give the group names
%   together with the axis; for instance, 'AB_1_1[0]' means the [0, 0, 1]
%   direction in the 'AB_1_1' plane. These are only used for labeling the
%   axes.

washold = ishold;
hold on;

angle_range = linspace(0, 2*pi, 100);

% draw circles for orientation
circle_radii = [1/4, 1/2, 3/4, 1];
for i = 1:length(circle_radii)
    crt_rad = circle_radii(i);
    plot(crt_rad*cos(angle_range), crt_rad*sin(angle_range), ':', 'color', [0.4 0.4 0.4]);
end

% draw the main axes
plot([-1 1], [0 0], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
plot([0 0], [-1 1], ':', 'color', [1 0.6 0.6], 'linewidth', 1);

% label the axes
xlabel(groupdir1);
ylabel(groupdir2);

if ~washold
    hold off;
end

end