function [Esbc] = set_thermalSBC(para)
% set the surface boundary condition for the thermal model
% 2019-2-15

global hS M

Cp = para.Cp;
Tref = para.Tref;

%% Storglaciaeren example
%
% T0 = 273.15;
% Tma =  -6.0; % degC, mean annual air temperature at Tarfala
% zcts = 1300; % m a.s.l.; altitude where CTS is at the surface, projected to topg
% slope = 100; % m; range around which surface temp transition happens
% Tsbc = T0 + Tma * (zcts + slope - hS) / (2.0 * slope);
% Tsbc(hS<zcts-slope) = T0 + Tma;
% Tsbc(hS>zcts+slope) = T0;

%% temperate ice experiment
T0 = 273.15 - 7.0; % 冰川末端年均气温
z0 = hS(end);

zfirn = median(hS); %; % 替换为每年的雪线高度
Tsbc = zeros(1,M);
for i=1:M
    if hS(i)<=zfirn
        Tsbc(i) = T0 - 6.5e-3*(hS(i)-z0); % 雪线高度以下，气温递减率
    else
        Tsbc(i) = 273.15 - 1.9; %273.15 - 2.5; % 雪线高度以上，假定都是-1.5℃
    end
end

%% Enthalpy
Esbc = Cp*(Tsbc - Tref);
end
