global xi N H zeta hS hB

figure

xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

contourf(xx/1000, yy, u, 30, 'LineStyle', 'none')
hold on
plot(xi/1000, hB, 'k', 'linewidth', 1);
plot(xi/1000, hS, 'k', 'linewidth', 1);
colormap('jet')

axis off

set(gcf,'position',[767.4000  341.8000  280.8000  168.8000])
