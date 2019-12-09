figure
global xi N H zeta

subplot(2,1,1)
xx = ones(N,1)*xi;
yy = ones(N,1)*hB(1:end-2) + zeta*H(1:end-2);

contourf(xx/1000, yy, u(:,1:end-1), 50, 'LineStyle', 'none')
hold on
plot(xi/1000, hB(1:end-2), 'k')
plot(xi/1000, hS(1:end-2), 'k')

c = colorbar;
set(c, 'Location', 'EastOutside', 'position', [0.91 0.62 0.02 0.25], 'fontsize', 7)
colormap('Jet')

ylabel(c, 'Velocity (m a^{-1})', 'fontsize', 7)

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')

subplot(2,1,2)
plot(xi/1000, u(end,1:end-1), '-')

xlabel('Horizontal distance (km)')
ylabel('Velocity (m a^{-1})')

title('Surface velocity')

set(gcf, 'Position', [423 220 553 433])