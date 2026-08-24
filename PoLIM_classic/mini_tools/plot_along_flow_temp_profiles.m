
addpath('./tools_plot')
global hS hB H xi zeta dzeta

icegrid = [9, 10];
textlabel = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'};
figure

%%
for i=1:12
    ith = icegrid(i);
    
    subplot(4,6,i)
    plot_temp_profile(H(ith), T(:,ith));
    text(0.1,0.1,textlabel{i},'unit','normalized','fontsize',12)
end

%%

subplot(4,6,13:18)

yyaxis left
plot(xi/1e3, hS, 'k', 'linewidth', 1);
hold on
plot(xi/1e3, hB, 'k', 'linewidth', 1);
for i=1:12
    ith = icegrid(i);
    plot([xi(ith)/1e3,xi(ith)/1e3], [hS(ith),hS(ith)-H(ith)], 'r-', 'linewidth', 3')
    text(xi(ith)/1e3+0.06, hS(ith)-H(ith)/2, textlabel{i})
end
xlabel('Distance (m)')
ylabel('Elevation (m)')
xlim([0, xi(end)/1e3+0.1])
ylim([hS(end)-30, hS(1)+30])

yyaxis right
Tsbc = Esbc/p.Cp + p.Tref;
plot(xi/1e3, Tsbc-273.15, 'linewidth', 1);
ylabel('T_{surf} (^\circC)')

xlim([0, xi(end)/1e3+0.1])

%%
ax = subplot(4,6,19:24);
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

hold on
contourf(xx/1000, yy, T-273.15, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
    plot(xi(index2)/1000, hB(index2)+CTS(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

cb = colorbar;
colormap(ax, 'jet')
pos_cb = get(cb,'Position');
set(cb, 'location','North', 'position',[pos_cb(1)-0.2,pos_cb(2)+0.07,0.27,0.02], 'fontsize', 10)
ylabel(cb, 'T (^\circC)', 'fontsize', 10)


xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
xlim([0, xi(end)/1e3+0.1])
ylim([hS(end)-30, hS(1)+30])
box on


set(gcf, 'position', 1e3*[0.2362,0.0674,1.0504,0.7016])
