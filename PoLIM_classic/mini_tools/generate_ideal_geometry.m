% Parameters
L_total = 10000;   % total length in meters
L_ice   = 8000;    % glacier-covered length in meters
H_max   = 200;     % maximum ice thickness in meters

% Create x grid
x = linspace(0, L_total, 51);  % in meters

% Bedrock: linearly from 300 m at head to 0 m at end
bedrock = 300 * (1 - x / L_total);

% Initialize surface as bedrock
surface = bedrock;

% Glacier mask
mask_ice = x <= L_ice;
x_ice = x(mask_ice);
x_norm = x_ice / L_ice;  % normalized [0,1]

% Parabolic ice thickness: 0 at both ends, max at center
thickness = H_max * 4 .* x_norm .* (1 - x_norm);

% Surface = bedrock + ice thickness in glacier region
surface(mask_ice) = bedrock(mask_ice) + thickness;

% Plot
figure('Position', [100, 100, 800, 400])
hold on;
fill_between = @(x, y1, y2, colorval) fill([x, fliplr(x)], [y1, fliplr(y2)], colorval, 'EdgeColor', 'none');
fill_between(x/1000, bedrock, surface, [0.8 0.92 1])  % ice fill
plot(x/1000, bedrock, 'Color', [0.55 0.27 0.07], 'LineWidth', 2)  % bedrock
plot(x/1000, surface, 'Color', [0 0.5 1], 'LineWidth', 2)         % surface
xline(L_ice/1000, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Glacier Terminus (8 km)', 'LabelVerticalAlignment', 'bottom')

xlabel('Distance along flowline (km)')
ylabel('Elevation (m)')
title('Idealized Glacier Geometry with Zero Thickness at Head and Terminus')
legend('Ice', 'Bedrock', 'Surface', 'Location', 'southwest')
grid on;
axis tight;
