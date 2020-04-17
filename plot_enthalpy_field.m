global xi N zeta dzeta H hS hB

figure
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

index1 = find(CTS>0 & CTS<N);
index2 = [index1(1)-1, index1, index1(end)+1];
% index2 = index1;

subplot(3,1,1)
hold on
contourf(xx/1000, yy, E, 20, 'LineStyle', 'none')
plot(xi/1000, hB, 'k')
plot(xi/1000, hS, 'k')
plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w--', 'linewidth',1.5)
hold off

c = colorbar;
colormap('Jet')

ylabel(c, 'Enthalpy (J kg^{-1})')

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')
box on

subplot(3,1,2)
hold on
contourf(xx/1000, yy, T-273.15, 20, 'LineStyle', 'none')
plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w--', 'linewidth',1.5)
plot(xi/1000, hB, 'k')
plot(xi/1000, hS, 'k')
hold off

c = colorbar;
colormap('Jet')

ylabel(c, 'Temperature (^\circC)')

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')
box on

subplot(3,1,3)
hold on
contourf(xx/1000, yy, omega*100, 20, 'LineStyle', 'none')
plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w--', 'linewidth',1.5)
plot(xi/1000, hB, 'k')
plot(xi/1000, hS, 'k')
hold off

c = colorbar;
%colormap('gray') % parula, jet, gray, winter
colormap(brewermap([],'*RdBu')) % *RdYlBu *RdBu

ylabel(c, 'Moisture (%)')

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')
box on

set(gcf, 'Position', [197 61 445 595])
