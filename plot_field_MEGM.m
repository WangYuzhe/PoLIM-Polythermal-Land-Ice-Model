% Date: 2018-2-2
% plot simulation results using MEGM scheme

global N xi zeta H hB hS SPY dzeta

SPY = 31556926;

% construct
xx = ones(N,1)*xi;
yy = zeta*H + ones(N,1)*hB;
index = find(CTS>0);

figure
subplot(3,2,1)
hold on
contourf(xx/1000, yy, E, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 1)
plot(xi/1000, hS, 'k', 'linewidth', 1)
plot(xi(index)/1000, H(index).*(CTS(index)-1)*dzeta + hB(index), 'w-', 'linewidth', 1)

xlabel('Distance (km)')
ylabel('Elevation (m asl)')
title('Enthalpy distribution')

c = colorbar;
title(c, 'J m^{-3}')
colormap('jet')
box on

subplot(3,2,2)
hold on
contourf(xx/1000, yy, T-273.15, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 1)
plot(xi/1000, hS, 'k', 'linewidth', 1)
plot(xi(index)/1000, H(index).*(CTS(index)-1)*dzeta + hB(index), 'w-', 'linewidth', 1)

xlabel('Distance (km)')
ylabel('Elevation (m asl)')
title('Temperature distribution')

c = colorbar;
title(c, '^\circC')
colormap('jet')
box on

subplot(3,2,3)
hold on
contourf(xx/1000, yy, omega*100, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 1)
plot(xi/1000, hS, 'k', 'linewidth', 1)
plot(xi(index)/1000, H(index).*(CTS(index)-1)*dzeta + hB(index), 'w-', 'linewidth', 1)

xlabel('Distance (km)')
ylabel('Elevation (m asl)')
title('Moisture distribution')

c = colorbar;
title(c, '%')
colormap('jet')
box on

subplot(3,2,4)
hold on
contourf(xx/1000, yy, u, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 1)
plot(xi/1000, hS, 'k', 'linewidth', 1)
xlabel('Distance (km)')
ylabel('Elevation (m asl)')
title('Velocity distribution')

c = colorbar;
title(c, 'm a^{-1}')
colormap('jet')
box on

subplot(3,2,[5,6])
hold on
plot(xi/1000, temperateWaterFlux(1, :)*SPY*1000, 'k')

xlabel('Distance (km)')
ylabel('Water flux (mm yr^{-1})')

ylim([-40 0])
box on

set(gcf, 'position', [206 84 933 534])