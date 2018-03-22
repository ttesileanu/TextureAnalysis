function drawTernaryTriangle
% drawTernaryTriangle Draw guiding simplex and circles for one of the
% ternary texture planes.

washold = ishold;
hold on;

max_t2 = sqrt(3)/2;
angle_range = linspace(0, 2*pi, 100);

% draw circles for orientation (radius 1 and 1/2)
plot(cos(angle_range), sin(angle_range), ':', 'color', [0.4 0.4 0.4]);
plot(0.5*cos(angle_range), 0.5*sin(angle_range), ':', 'color', [0.4 0.4 0.4]);

% draw the probability triangle
plot([-1/2 1 -1/2 -1/2], [-max_t2 0 max_t2 -max_t2], 'color', [0.5 0.7 1]);

% draw the main axes
plot([0 1.5], [0 0], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
plot([0 -1/2*1.5], [0 1.5*max_t2], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
plot([0 -1/2*1.5], [0 -1.5*max_t2], ':', 'color', [1 0.6 0.6], 'linewidth', 1);

% label the corners
text(1.07, -0.1, '[0,1,0]', 'fontsize', 12, 'color', [0.7 0.7 0.7]);
text(-1.1,  max_t2+0.05, '[0,0,1]', 'fontsize', 12, 'color', [0.7 0.7 0.7]);
text(-1.1, -max_t2-0.01, '[1,0,0]', 'fontsize', 12, 'color', [0.7 0.7 0.7]);

if ~washold
    hold off;
end

end