% 2019-1-16 
function p = params_Kleiner_ExpB()

p.SPY    = 31556926;    % seconds in one year [s a-1]
p.rhoi   = 910;         % ice density [kg m-3]
p.rhow   = 1000;        % water density [kg m-3]
p.g      = 9.81;        % gravitational acceleration [m s-2]
p.kc     = 2.1;         % cold ice conductivity [W m-1 K-1]
p.Cp     = 2009;        % ice specific heat capacity [J kg-1 K-1]
p.Lw     = 3.35e5;      % latent heat of fusion [J kg-1]
p.Tref   = 223.15;      % Reference temperature [K] 273, 223.15;
p.AGlen  = 5.3e-24;     % rate factor [Pa-3 s-1]
p.Kc     = p.kc/(p.rhoi*p.Cp); % Thermal diffusivity of cold ice [m2 s-1]
p.Kt     = p.Kc*1e-5;      % Thermal diffusivity of temperate ice [m2 s-1]
p.k0 = 1e-12; % permeability factor (unconstrained, Hewitt&Schoof, 2017, Tab. 1) [m2]; range: 1e-12~5e-8
p.eta_w  = 1.8e-3;      % viscosity of water [Pa s]
p.alpha_hewitt = 2; % exponent in compaction pressure model (Hewitt&Schoof, 2017); unconstrained, range: 2~3
p.betaCC = 0;           % Clausius-Clapeyron constant [K Pa-1]
p.Qgeo   = 0;
p.n      = 3;

p.type_enth_solver = 'SEGM';
p.type_valley = 'rect';
p.Hmin = 0.0;
p.is_greve_drain = 0;
p.is_auto_enth_BBC = 1;
p.type_enth_BBC = 'TEMP_base';
p.layers = 201;
p.is_enth_trans = 1;
end