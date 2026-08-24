% Date: 2017-7-23
global kc rho g Cp Lw Qgeo SPY betaCC Tref rhow AGlen Kc Kt type_geometry

SPY    = 31556926;    % seconds in one year [s a-1]
rho    = 910;         % ice density [kg m-3]
rhow   = 1000;        % water density [kg m-3]
g      = 9.81;        % gravitational acceleration [m s-2]
kc     = 2.1;         % cold ice conductivity [W m-1 K-1]
Cp     = 2009;        % ice specific heat capacity [J kg-1 K-1]
Lw     = 3.35e5;      % latent heat of fusion [J kg-1]
Qgeo   = 0;           % geothermal heat flux [W m-2] 
betaCC = 0;           % Clausius-Clapeyron constant [K Pa-1]
Tref   = 223.15;      % Reference temperature [K]
AGlen  = 5.3e-24;     % rate factor [Pa-3 s-1]
Kc     = kc/(rho*Cp); % Thermal diffusivity of cold ice [m2 s-1]
Kt     = 1.1e-11;      % Thermal diffusivity of temperate ice [m2 s-1][m2 s-1]

type_geometry = 4;
