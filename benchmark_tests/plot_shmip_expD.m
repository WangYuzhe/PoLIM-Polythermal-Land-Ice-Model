% time
SPD = 24*3600;
SPY = 365*SPD;

dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

timeSpan = arrayTime/SPD;

% plot
figure
subplot(5,1,1)
load result_D1.mat
plot_benchmark.plot_shmip_D(timeSpan, At_N)

subplot(5,1,2)
load result_D2.mat
plot_benchmark.plot_shmip_D(timeSpan, At_N)

subplot(5,1,3)
load result_D3.mat
plot_benchmark.plot_shmip_D(timeSpan, At_N)

subplot(5,1,4)
load result_D4.mat
plot_benchmark.plot_shmip_D(timeSpan, At_N)

subplot(5,1,5)
load result_D5.mat
plot_benchmark.plot_shmip_D(timeSpan, At_N)

axis tight
set(gcf, 'position', [521.0696   78.6348  521.7391  720.4174])