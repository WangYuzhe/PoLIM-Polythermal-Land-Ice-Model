function [ax] = plot_temp_profile(thk,iceT)
%PLOT_TEMP_PROFILE Summary of this function goes here
%   Detailed explanation goes here
global zeta

depth = thk*zeta;

iceT = iceT - 273.15;
iceT = flipud(iceT);

% ax = plot(iceT, depth, '-o', 'linewidth', 1.5, 'markersize', 3);
ax = plot(iceT, depth, '-', 'linewidth', 1.5);

xlabel('Temperature (^\circC)')
ylabel('Depth (m)')

set(gca, 'YDir', 'reverse')

end

