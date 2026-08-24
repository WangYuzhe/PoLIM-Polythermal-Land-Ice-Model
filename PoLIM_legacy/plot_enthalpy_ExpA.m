load enthA_analy_result.mat

figure
subplot(4,1,1)
plot(arrayTime/1000, At_basalT(:,1)-273.15, 'k', 'linewidth', 1)
xlim([0 310])
ylim([-30 5])
ylabel('T_b (^\circC)')
grid on
box on

subplot(4,1,2)
hold on
plot(basalMelt(2,:)/1000, basalMelt(1,:), 'r', 'linewidth', 2)
plot(arrayTime/1000, At_basalMeltRate(:,1)*1000, 'k', 'linewidth', 1)
xlim([0 310])
ylim([-4 4])
ylabel('a_b (mm a^{-1} w.e.)')
grid on
box on

subplot(4,1,3)
plot(arrayTime/1000, At_Hw(:,1), 'k', 'linewidth', 1)
xlim([0 310])
ylim([-10 160])
ylabel('H_w (m)')
xlabel('Time (ka)')
grid on
box on

subplot(4,1,4)
plot(arrayTime/1000, At_isTemperate(:,1), 'k', 'linewidth', 1)
xlim([0 310])
ylim([-0.2 1.2])
ylabel('Basal thermal state')
xlabel('Time (ka)')
grid on
box on

set(gcf, 'position', [246 114 831 536])