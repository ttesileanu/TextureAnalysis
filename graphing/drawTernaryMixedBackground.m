function drawTernaryMixedBackground(groupDir1, groupDir2)
% drawTernaryMixedBackground Draw guiding axes through origin and guiding
% circles for pairs of ternary texture groups.
%   drawTernaryMixedBackground(groupDir1, groupDir2) draws guiding axes
%   going through the origin and guiding circles for a plot representing a
%   mix of two ternary texture groups. The arguments give the group names
%   together with the axis; for instance, a typical notation would write
%   'AB_1_1[0]' for the [0, 0, 1] direction in the 'AB_1_1' plane. However,
%   these are only used for labeling the axes, so any other notation can be
%   used.

% we will overlay several graphical objects
wasHold = ishold;
hold on;

% we use this for drawing the circles
angleRange = linspace(0, 2*pi, 100);

% draw circles for orientation
circleRadii = [1/4, 1/2, 3/4, 1];
for i = 1:length(circleRadii)
    radiuds = circleRadii(i);
    plot(radiuds*cos(angleRange), radiuds*sin(angleRange), ':', 'color', [0.4 0.4 0.4]);
end

% draw the main axes
plot([-1 1], [0 0], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
plot([0 0], [-1 1], ':', 'color', [1 0.6 0.6], 'linewidth', 1);

% label the axes
xlabel(groupDir1);
ylabel(groupDir2);

% revert hold state
if ~wasHold
    hold off;
end

end