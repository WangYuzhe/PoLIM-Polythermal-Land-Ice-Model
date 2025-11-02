% Date: 2018-2-2
% plot simulation results using SEGM scheme (without Greve drainage function)

global N xi zeta H hB hS SPY dzeta

SPY = 31556926;

index1 = find(CTS>0 & CTS<N);

if ~isempty(index1)
    if index1(1)==1
        if index1(end)==M
            index2 = index1; % 1:M
        else
            index2 = [index1, index1(end)+1];
        end
    else
        if index1(end)==M
            index2 = [index1(1)-1, index1];
        else
            index2 = [index1(1)-1, index1, index1(end)+1];
        end
    end
end


% construct
xx = ones(N,1)*xi;
yy = zeta*H + ones(N,1)*hB;

figure
subplot(3,2,1)
hold on
contourf(xx/1000, yy, E, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 1)
plot(xi/1000, hS, 'k', 'linewidth', 1)
if ~isempty(index1)
    plot(xi(index2)/1000, H(index2).*(CTS(index2)-1)*dzeta + hB(index2), 'w-', 'linewidth', 1.5)
end

xlabel('Distance (km)')
ylabel('Elevation (m asl)')
title('Enthalpy distribution')

c = colorbar;
title(c, 'J kg^{-1}')
colormap('jet')
box on

subplot(3,2,2)
hold on
contourf(xx/1000, yy, T-273.15, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 1)
plot(xi/1000, hS, 'k', 'linewidth', 1)

if ~isempty(index1)
    plot(xi(index2)/1000, H(index2).*(CTS(index2)-1)*dzeta + hB(index2), 'w-', 'linewidth', 1.5)
end

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

if ~isempty(index1)
    plot(xi(index2)/1000, H(index2).*(CTS(index2)-1)*dzeta + hB(index2), 'w-', 'linewidth', 1.5)
end

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

if strcmpi(p.type_enth,'SEGM')
    if p.is_greve_drain
        plot(xi/1000, -qw_TEMP_drain_greve*SPY*1000, 'k')
        min_qw_TEMP_drain_greve = min(-qw_TEMP_drain_greve*SPY*1000);
        ylim([floor(min_qw_TEMP_drain_greve)-5, 0])
    else
        plot(xi/1000, -qw_TEMP_diffu(2,:)*SPY*1000, 'k')
    end
else
    plot(xi/1000, qw_TEMP(1,:)*SPY*1000, 'k-')
    hold on
    plot(xi/1000, qw_TEMP_darcy(1,:)*SPY*1000, 'r-')
    ylim([-45 0])
    legend('Total flux', 'Darcy flux', 'Location', 'SouthWest')
end

xlabel('Distance (km)')
ylabel('Water flux (mm yr^{-1})')

box on

set(gcf, 'position', [206 84 933 534])