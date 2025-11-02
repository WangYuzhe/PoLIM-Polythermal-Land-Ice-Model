function [strainHeat] = calc_strainHeat()
% Date: 2025/10/31

strain_heat_s = 4*visc_s.*de2_s; % [Pa a-1]
strain_heat = staggerX2main(strain_heat_s); % [Pa a-1]

end