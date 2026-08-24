global xi N H zeta hS hB

%% plot1
ax = subplot(2,1,1);
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

contourf(xx/1000, yy, u, 30, 'LineStyle', 'none')
hold on
plot(xi/1000, hB, 'k', 'linewidth', 1);
plot(xi/1000, hS, 'k', 'linewidth', 1);

cb = colorbar(ax);
pos_cb = get(cb,'Position');
set(cb, 'location','North', 'position',[pos_cb(1)-0.2,pos_cb(2)+0.23,0.27,0.02], 'fontsize', 10)
% set(cb, 'location','North', 'position',[pos_cb(1)-0.1,pos_cb(2)+0.3,0.15,0.02], 'fontsize', 10)
ylabel(cb, 'Velocity (m a^{-1})', 'fontsize', 10)

colormap('jet')

xlabel('Distance from headwall (km)', 'fontsize', 12)
ylabel('Elevation (m a.s.l.)', 'fontsize', 12)

%% plot2
subplot(2,1,2)
HL = zeros(1,2);
HL(1) = plot(xi/1000, u(end, :), 'k-', 'linewidth', 1.5);
hold on
HL(2) = plot(xi/1000, u(1,:), 'k--', 'linewidth', 1.5);

xlabel('Distance from headwall (km)', 'fontsize', 12)
ylabel('Velocity (m a^{-1})', 'fontsize', 12)

hlgd = legend(HL, 'Surface velocity', 'Sliding velocity');
set(hlgd, 'fontsize', 9, 'location', 'NorthEast', 'color','none');

set(gcf, 'Position', [507,110,469,489])
