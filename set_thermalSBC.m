function [Esbc] = set_thermalSBC()
% set the surface boundary condition for the thermal model
% 2019-2-15

global Cp Tref hS hB M

% Storglaciaeren example
T0 = 273.15;
Tma =  -6.0; % degC, mean annual air temperature at Tarfala
zcts = 1300; % m a.s.l.; altitude where CTS is at the surface, projected to topg
slope = 100; % m; range around which surface temp transition happens
Tsbc = T0 + Tma * (zcts + slope - hS) / (2.0 * slope);
Tsbc(hS<zcts-slope) = T0 + Tma;
Tsbc(hS>zcts+slope) = T0;
Esbc = Cp*(Tsbc - Tref);
end
