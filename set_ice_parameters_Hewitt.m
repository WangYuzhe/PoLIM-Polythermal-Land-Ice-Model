% 2019-1-16 
global kc rho g Cp Lw SPY Tref rhow AGlen Kc Kt eta_w k0 betaCC Qgeo de0 n type_geometry

SPY    = 31556926;    % seconds in one year [s a-1]
rho    = 916;         % ice density [kg m-3]
rhow   = 1000;        % water density [kg m-3]
g      = 9.81;        % gravitational acceleration [m s-2]
kc     = 2.1;         % cold ice conductivity [W m-1 K-1]
Cp     = 2009;        % ice specific heat capacity [J kg-1 K-1]
Lw     = 3.34e5;      % latent heat of fusion [J kg-1]
Tref   = 223.15;      % Reference temperature [K]
AGlen  = 2.4e-24;     % rate factor [Pa-3 s-1]
Kc     = kc/(rho*Cp); % Thermal diffusivity of cold ice [m2 s-1]
Kt     = 1.1e-8;      % Thermal diffusivity of temperate ice [m2 s-1]
eta_w  = 1.8e-3;      % viscosity of water [Pa s]
k0     = 1e-12;       % permeability factor [m2]
betaCC = 0;           % Clausius-Clapeyron constant [K Pa-1]
Qgeo   = 0;
de0    = 1e-30; % [a-2]; small number in case of singularity
n      = 3;

type_geometry = 4;