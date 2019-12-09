clear, clc

load N_E1.mat
load N_E2.mat
load N_E3.mat
load N_E4.mat
load N_E5.mat

dx = 60;
x = 0:dx:6000;
M = length(x);

y1 = 1*ones(1,M);
y2 = 2*ones(1,M);
y3 = 3*ones(1,M);
y4 = 4*ones(1,M);
y5 = 5*ones(1,M);

figure
hold on
plot3(x/1000, y1, fliplr(N_E1/1e6), 'k', 'linewidth', 2)
plot3(x/1000, y2, fliplr(N_E2/1e6), 'k', 'linewidth', 2)
plot3(x/1000, y3, fliplr(N_E3/1e6), 'k', 'linewidth', 2)
plot3(x/1000, y4, fliplr(N_E4/1e6), 'k', 'linewidth', 2)
plot3(x/1000, y5, fliplr(N_E5/1e6), 'k', 'linewidth', 2)
grid on

xlabel('Distance (km)')
zlabel('N (MPa)')
zlim([0, 4])
set(gca, 'YDir', 'reverse')

set(gcf, 'position', [460 101 744 553])

view(-120, 22)