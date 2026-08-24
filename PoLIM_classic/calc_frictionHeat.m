function [frictionHeat] = calc_frictionHeat(u_s, u, visc_s, para)
% Date: 2023/12/2

% calculate the friction heat
[tauxz, ~] = calc_tauxz(u_s, visc_s);
frictionHeat = tauxz.*u(1,:)/para.SPY; % [W m-2]

end