figure
global xi N H zeta hS hB

subplot(2,1,1)
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

contourf(xx/1000, yy, u, 50, 'LineStyle', 'none')
hold on
plot(xi/1000, hB, 'k')
plot(xi/1000, hS, 'k')

c = colorbar;
% set(c, 'Location', 'North', 'position', [0.62 0.785 0.25 0.02], 'fontsize', 9)
set(c, 'Location', 'EastOutside', 'position', [0.91 0.62 0.02 0.25], 'fontsize', 9)

colormap('Jet')

ylabel(c, 'Velocity (m a^{-1})')

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')

subplot(2,1,2)
plot(xi/1000, u(end, :), 'k-')

xlabel('Horizontal distance (km)')
ylabel('Velocity (m a^{-1})')

% ylim([0 max(max(u)) + 20])

title('Surface velocity')

set(gcf, 'Position', [423 145 689 508])