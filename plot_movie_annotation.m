% Plot transient velocity results as a movie with max/min annotations

global xi N hS hB zeta

xx = ones(N,1) * xi / 1e3;  % Grid in km

zmax = ceil(max(hS)) + 100;
zmin = floor(min(hB)) - 100;

% Prepare figure
figure(1)
set(gcf, 'Position', [200, 200, 800, 400])

for i = 1:numTimeStep
    clf;

    % Elevation grid at time i
    yy = ones(N,1) * hB + zeta * At_H(i,:);
    u_yr = At_u(:,:,i);              % Velocity field at time i
    umax = max(u_yr(:));
    umin = min(u_yr(:));

    % Indices of max and min velocity
    [max_row, max_col] = wyz_max_index(u_yr);
    [min_row, min_col] = wyz_min_index(u_yr);

    % Physical locations of max/min
    x_umax = xi(max_col)/1e3;
    z_umax = hB(max_col) + (max_row - 1)/(N - 1) * H(max_col);

    x_umin = xi(min_col)/1e3;
    z_umin = hB(min_col) + (min_row - 1)/(N - 1) * H(min_col);

    % Plotting
    hold on
    contourf(xx, yy, u_yr, 30, 'LineStyle', 'none')
    plot(xi / 1e3, hB, 'k', 'LineWidth', 1.5)
    plot(xi / 1e3, At_hS(i,:), 'k', 'LineWidth', 1.5)
    plot(x_umax, z_umax, 'ro', 'MarkerSize', 5, 'LineWidth', 1.5)
    plot(x_umin, z_umin, 'k^', 'MarkerSize', 5, 'LineWidth', 1.5)

    % Colorbar and labels
    ch = colorbar('EastOutside');
    set(get(ch, 'Title'), 'String', '[m/a]', 'FontSize', 14);
    colormap('jet')

    % Add text annotations for max and min values
    text(x_umax + 0.05, z_umax+40, sprintf('max: %.2f', umax), ...
        'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center')

    text(x_umin + 0.05, z_umin-30, sprintf('min: %.2f', umin), ...
        'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center')

    % Axes and title
    ylim([zmin, zmax])
    xlim([xi(1)/1e3-0.1, xi(end)/1e3+0.5])

    xlabel('x (km)', 'FontSize', 14)
    ylabel('z (m)', 'FontSize', 14)
    title(['Model Year: ', num2str(arrayTime(i), '%10.0f')], 'FontSize', 16)
    set(gca, 'FontSize', 14)
    box on

    % Capture frame
    F(i) = getframe(gcf);
end

% % Write to video
% outfile = './outputs/ideal_glac_transient_thk_evol_30yrs.avi';
% videoFile = VideoWriter(outfile);
% videoFile.FrameRate = 4;
% open(videoFile);
% writeVideo(videoFile, F);
% close(videoFile);
