load Rechg_D3.mat
load result_D3.mat

% PARAMETERS
SPD = 24*3600; % [s d-1]
SPY = 365*SPD; % [s yr-1]

% TIME SETTING
dt_hydro = 1/24*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
timeSpan = 0:dt_hydro:(endTime-dt_hydro);
numTimeStep = length(timeSpan);

figure
hold on
plot(timeSpan/SPD, At_Qw(:,1), 'r--', 'linewidth', 1.5);
plot(timeSpan/SPD, At_Qw(:,2), 'c--', 'linewidth', 1.5);
plot(timeSpan/SPD, At_Qw(:,3), 'b--', 'linewidth', 1.5);

plot(timeSpan/SPD, bandRechg(:,1), 'r-', 'linewidth', 1.5);
plot(timeSpan/SPD, bandRechg(:,2), 'c-', 'linewidth', 1.5);
plot(timeSpan/SPD, bandRechg(:,3), 'b-', 'linewidth', 1.5);

% plot(timeSpan/SPD, totalRechg(:,1), 'k-', 'linewidth', 1.5);

hold off
box on
xlim([90, 365])
set(gca, 'xtick', [120 181 243 304 365])
set(gca, 'xticklabel', [4 6 8 10 12])

xlabel('Time (month)')
ylabel('Discharge (m^3 s^{-1})')
