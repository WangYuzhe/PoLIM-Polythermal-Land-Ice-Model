load enthB_analy_result.mat

figure
subplot(1,3,1)
hold on
plot(enthB_analy_E/1000, enthB_analy_z, 'r', 'linewidth', 2)
plot(At_E(:,1,end)/1000, zeta, '-k', 'linewidth', 1)
% plot(At_E(jcts,:,end)/10^3, zeta(jcts), 'ob')
xlim([94 108])
set(gca, 'xtick', [96 100 104 108])
set(gca, 'ytick', 0:0.1:1)

xlabel('E (\times 10^3 J kg^{-1})')
ylabel('\zeta')
grid on

subplot(1,3,2)
hold on
plot(enthB_analy_T-273.15, enthB_analy_z, 'r', 'linewidth', 2)
plot(At_T(:,1,end)-273.15, zeta, '-k', 'linewidth', 1)
set(gca, 'xtick', -3:0.5:0.5)
set(gca, 'ytick', 0:0.1:1)
set(gca, 'YTickLabel', [])

xlim([-3 0.5])
set(gca, 'xtick', [-3 -2 -1 0])
xlabel('T (^\circC)')
grid on

subplot(1,3,3)
hold on
plot(enthB_analy_omega*100, enthB_analy_z, 'r', 'linewidth', 2)
plot(At_omega(:,1,end)*100, zeta, '-k', 'linewidth', 1)
set(gca, 'xtick', 0:0.5:3)
set(gca, 'ytick', 0:0.1:1)
set(gca, 'YTickLabel', [])

xlim([-0.5 3])
xlabel('\omega (%)')
grid on

set(gcf, 'position', [239 153 1005 412])