global xi N H zeta hS hB

figure
subplot(2,1,1)
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

contourf(xx/1000, yy, u, 50, 'LineStyle', 'none')
hold on
plot(xi/1000, hB, 'k')
plot(xi/1000, hS, 'k')

vcb = colorbar;
set(vcb, 'location', 'North', 'position', [0.65 0.8 0.25 0.02], 'fontsize', 8)
colormap('Jet')

ylabel(vcb, 'Velocity (m a^{-1})', 'fontsize', 8)

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')

subplot(2,1,2)
HL = zeros(1,2);
HL(1) = plot(xi/1000, u(end, :), 'b-');
hold on
HL(2) = plot(xi/1000, u(1,:), 'r--');

xlabel('Horizontal distance (km)')
ylabel('Velocity (m a^{-1})')

hlgd = legend(HL, 'surface velocity', 'sliding velocity');
set(hlgd, 'fontsize', 8, 'location', 'NorthWest');

ylim([0 max(max(u)) + 20])
title('Surface velocity')

set(gcf, 'Position', [581 147 455 530])