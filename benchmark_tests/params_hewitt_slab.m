% 2019-1-16 
function p = params_hewitt_slab()

p.SPY    = 31556926;    % seconds in one year [s a-1]
p.rhoi   = 916;         % ice density [kg m-3]
p.rhow   = 1000;        % water density [kg m-3]
p.g      = 9.81;        % gravitational acceleration [m s-2]
p.kc     = 2.1;         % cold ice conductivity [W m-1 K-1]
p.Cp     = 2009;        % ice specific heat capacity [J kg-1 K-1]
p.Lw     = 3.34e5;      % latent heat of fusion [J kg-1]
p.Tref   = 273.00;      % Reference temperature [K] 273, 223.15;
p.AGlen  = 2.4e-24;     % rate factor [Pa-3 s-1]
p.Kc     = p.kc/(p.rhoi*p.Cp); % Thermal diffusivity of cold ice [m2 s-1]
p.Kt     = 1.1e-8;      % Thermal diffusivity of temperate ice [m2 s-1]
p.eta_w  = 1.8e-3;      % viscosity of water [Pa s]
p.k0     = 1e-12;       % permeability factor [m2]
p.betaCC = 0;           % Clausius-Clapeyron constant [K Pa-1]
p.Qgeo   = 0;
p.de0    = 1e-30; % [a-2]; small number in case of singularity
p.n      = 3;

p.type_enth_solver = 'MEGM';
p.type_valley = 'rect';
p.Hmin = 1.0;

p.is_auto_enth_BBC = 0;
p.type_enth_BBC = 2;
p.has_Greve_drainage = 1;
p.layers = 201;
p.alpha_hewitt = 2; % exponent in compaction pressure model (Hewitt&Schoof, 2017); unconstrained, range: 2~3

p.is_enth_trans = 1;
% true (1): transient enthalpy solver
% false (0): steady-state enthalpy solver

end