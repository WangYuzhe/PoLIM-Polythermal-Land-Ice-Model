%==============================FIGURE==============================
figure
subplot(2,1,1)
hold on
FH = zeros(3,1);
FH(1) = plot(arrayTime/SPD, mean(At_N(:, Hindex1:Hindex2),2)/1e6, 'r', 'linewidth', 1.5);
FH(2) = plot(arrayTime/SPD, mean(At_N(:, Mindex1:Mindex2),2)/1e6, 'c', 'linewidth', 1.5);
FH(3) = plot(arrayTime/SPD, mean(At_N(:, Lindex1:Lindex2),2)/1e6, 'b', 'linewidth', 1.5);
hold off
box on

xlim([90, 365])
switch type_glacier_geometry
    case 2 % ice sheet
        ylim([-4.5, 10])
    case 3 % valley glacier
        ylim([-0.5, 6])
end
set(gca, 'xtick', [120 181 243 304 365])
set(gca, 'xticklabel', [4 6 8 10 12])
set(gca, 'ytick', [-4 0 4 8])
set(gca, 'yticklabel', [-4 0 4 8])

legend(FH, 'upper band', 'middle band', 'lower band', 'location', 'SouthEast')
xlabel('Time (month)')
ylabel('Effective pressure (MPa)')

subplot(2,1,2)
hold on
plot(arrayTime/SPD, At_Qw(:,1), 'r--', 'linewidth', 1.5);
plot(arrayTime/SPD, At_Qw(:,2), 'c--', 'linewidth', 1.5);
plot(arrayTime/SPD, At_Qw(:,3), 'b--', 'linewidth', 1.5);
hold off
box on
xlim([90, 365])
set(gca, 'xtick', [120 181 243 304 365])
set(gca, 'xticklabel', [4 6 8 10 12])

xlabel('Time (month)')
ylabel('Discharge (m^3 s^{-1})')
%==============================FIGURE==============================