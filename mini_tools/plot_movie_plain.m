% Plot transient velocity results as a movie with max/min annotations

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);

addpath(repo_root)
addpath(script_dir)

%% unpack
global xi zeta hB

%% initialization
glac_name = "Arolla";
outdir = fullfile(repo_root, "mov_results");
if ~isfolder(outdir)
    mkdir(outdir)
end

outfile = fullfile(outdir, "mov_" + glac_name);

[N, M, nt] = size(At_u);

xx = ones(N,1) * xi / 1e3;  % Grid in km

zmax = ceil(max(At_hS(1,:))) + 100;
zmin = floor(min(hB)) - 100;

umax = ceil(max(At_u(:)));

% --- Preallocate video ---
videoFile = VideoWriter(outfile);
videoFile.FrameRate = 4;
open(videoFile);

figure(1)

for i = 1:nt
    clf;

    % Elevation grid at time i
    yy = ones(N,1) * hB + zeta * At_H(i,:);
    u_yr = At_u(:,:,i); % Velocity field at time i

    % Plotting
    hold on
    contourf(xx, yy, u_yr, 30, 'LineStyle', 'none')
    plot(xi / 1e3, hB, 'k', 'LineWidth', 1)
    plot(xi / 1e3, At_hS(i,:), 'k', 'LineWidth', 1)

    % Colorbar and labels
    ch = colorbar('EastOutside');
    set(get(ch, 'Title'), 'String', '[m/a]', 'FontSize', 14);
    colormap('jet')
    % clim([0, umax]);

    % Axes and title
    ylim([zmin, zmax])
    xlim([xi(1)/1e3 - 0.1, xi(end)/1e3 + 0.2])

    xlabel('Distance from headwall (km)', 'fontsize', 14)
    ylabel('Elevation (m a.s.l.)', 'fontsize', 14)

    title(['Model Year: ', num2str(arrayTime(i), '%10.0f')], 'FontSize', 16)
    set(gca, 'FontSize', 14)
    set(gcf,'position',[200, 200, 800, 400])
    box on

    % --- Capture & write frame immediately (no large F array) ---
    frame = getframe(gcf);
    writeVideo(videoFile, frame);
end
hold off

close(videoFile);
