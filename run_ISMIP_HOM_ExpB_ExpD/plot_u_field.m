figure
global xi N H zeta

xx = ones(N,1)*xi;
yy = ones(N,1)*hB(1:end-2) + zeta*H(1:end-2);

figure
contourf(xx/1000, yy, u(:,1:end-1), 50, 'LineStyle', 'none')
hold on
plot(xi/1000, hB(1:end-2), 'k')
plot(xi/1000, hS(1:end-2), 'k')

c = colorbar;
colormap('Jet')

ylabel(c, 'Velocity (m a^{-1})')

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')

set(gcf, 'Position', [300 220 676 380])