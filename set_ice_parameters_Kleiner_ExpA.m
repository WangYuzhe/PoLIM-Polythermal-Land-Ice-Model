% Date: 2017-7-23
global kc rho g Cp Lw Qgeo SPY betaCC Tref rhow Kc Kt type_geometry

kc = 2.1;         % cold ice conductivity [W m-1 K-1]
rho = 910;        % ice density [kg m-3]
g = 9.81;         % gravitational acceleration [m s-2]
Cp = 2009;        % ice specific heat capacity [J kg-1 K-1]
Lw = 3.34e5;      % latent heat of fusion [J kg-1]
Qgeo = 0.042;        % geothermal heat flux [W m-2] 
SPY = 31556926;   % seconds in one year [s a-1]
betaCC = 7.9e-8;  % Clausius-Clapeyron constant [K Pa-1]
Tref = 223.15;    % Reference temperature [K]
rhow = 1000;      % water density [kg m-3]

Kc = kc/(rho*Cp); % [m2 s-1]
Kt = 1.1e-9; % [m2 s-1]

type_geometry = 3;
