% Date: 2017-9-11

function [tauxz, tauxz_s] = calc_tauxz(u_s, visc_s)
% tau_xz = eta / H * du/dzeta

global H_s dzeta

tauxz_s = visc_s(1,:)./H_s.*(-3*u_s(1,:)+4*u_s(2,:)-u_s(3,:))/(2*dzeta);

tauxz_s = tauxz_s/1000;

tauxz = [tauxz_s(:,1), (tauxz_s(:, 2:end-2)+tauxz_s(:, 3:end-1))/2, tauxz_s(:,end)];

end