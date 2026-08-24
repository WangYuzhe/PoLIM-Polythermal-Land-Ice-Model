% Date: 2023/11/7

figure
load('../geo_inputs/geo_arolla.mat')
f(1) = plot(geo.xi/1e3, geo.hS, 'b', 'linewidth', 1.5);
hold on
f(2) = plot(geo.xi/1e3, geo.hB, 'k-', 'linewidth', 1.5);

xtext = [0.50,0.58];
ytext = [0.4,0.4];

annotation('textarrow',xtext,ytext,'String','Thickness (H)')
idx = 30;
plot([geo.xi(idx)/1e3, geo.xi(idx)/1e3],[geo.hB(idx),geo.hS(idx)], 'r-', 'linewidth', 1.5);

xlabel('Distance from headwall (km)')
ylabel('Elevation (m asl.)')

legend(f, 'Surface (hS)', 'Bedrock (hB)')

