plot(xi/1e3, qw_TEMP_darcy(1,:)*SPY*1e3, 'k-')
hold on
plot(xi/1e3, -m_basal*1e3, 'r-')
xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Water flux (mm a^{-1})', 'FontSize', 10)
box on
grid on

legend('Darcy flux', 'Basal melt', 'Location', 'SouthWest')
