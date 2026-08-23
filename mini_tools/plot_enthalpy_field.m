global xi M N zeta dzeta H hS hB
SPY = 31536000;
  
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

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

%% Enthalpy field
%
ax(1) = subplot(2,2,1);
hold on
contourf(xx/1000, yy, E, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
    plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(1), 'jet');

ylabel(c, 'Enthalpy (J kg^{-1})')

ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
xlim([0, xi(end)/1e3])
box on
title('Enthalpy field')

%% Porosity field
%
ax(2) = subplot(2,2,2);
hold on
contourf(xx/1000, yy, omega*100, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(2), 'jet');

ylabel(c, 'Moisture (%)')

xlim([0, xi(end)/1e3])

box on

title('Moisture field')

%%
ax(3) = subplot(2,2,3);
hold on
contourf(xx/1000, yy, T-273.15, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
    plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(3), 'jet');

ylabel(c, 'Temperature (^\circC)')

xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
xlim([0, xi(end)/1e3])
box on
title('Temperature field')

%% Velocity
%
ax(4) = subplot(2,2,4);
hold on
contourf(xx/1000, yy, u, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(4), 'jet')

ylabel(c, 'u (m a^{-1})')

xlabel('Horizontal distance (km)', 'FontSize', 10)
xlim([0, xi(end)/1e3])

box on

title('Horizontal velocity field')

set(gcf, 'Position', [558,275,831,467])
