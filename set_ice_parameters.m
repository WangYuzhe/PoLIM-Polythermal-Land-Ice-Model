global g R SPY SPD rho rhow n de0 kc Cp Lw betaCC Qgeo Tref Kc Kt eta_w k0 epsilon iter_max 
global type_BBC type_valley isFlowband type_Arrhenius type_thermal_model type_geometry
global lambda_max m_max epr

% fixed physical constants
g = 9.81; % accel of gravity [m s-2]
R = 8.314; % universal gas constant [J mol K-1]; 
SPY = 31556926; % year is this many seconds (i.e. 365.2422 days)
SPD = 86400; % day is this many seconds
rho = 910; % density of ice [kg m-3]
rhow = 1000; % water density [kg m-3]
n = 3; % exponent in Glen's flow law []
de0 = 1e-30; % small number in case of singularity [yr-2]
kc = 2.1; % cold ice conductivity [W m-1 K-1]
Cp = 2009; % ice specific heat capacity [J kg-1 K-1]
Lw = 3.34e5; % latent heat of fusion [J kg-1]
betaCC = 7.9e-8; % Clausius-Clapeyron constant [K Pa-1]
Qgeo = 0.042; % geothermal heat flux [W m-2] 
Tref = 223.15; % Reference temperature [K]
eta_w = 1.8e-3; % viscosity of water [Pa s]
k0 = 1e-12; % permeability factor [m2]

Kc = kc/(rho*Cp); % [m2 s-1]
Kt = 1.1e-9; % [m2 s-1]

% adjustable parameters
epsilon = 0; % enhanced stress-free condition
iter_max = 50; % iterations

lambda_max = 4; % ref: 4
m_max = 0.3; % ref: 0.3
epr = 1; % N=Pi-Pw, epr=0, Pw=Pi; epr=1, Pw=0.

type_geometry = 6;
% 1: longitudinal profile of mountain glacier (default is Haut d'Arolla)
% 2: Kleiner Exp. A (parallel-sided slab)
% 3: Kleiner Exp. B (polythermal parallel-sided slab)
% 4: Hewitt-Schoof Ice Cap Experiment
% 5: SHMIP valley glacier

type_BBC = 2;
% 1: no-slip bed
% 2: Coulomb friction law
% 3: linear friction law

type_valley = 3;
% 1: trapezoid
% 2: Svessen: y=ax^b
% 3: rectangular

type_Arrhenius = 2;
% 1: Greve (2009)
% 2: Cuffey & Paterson (2010)

type_thermal_model = 1;
% 1: standard enthalpy gradient model (SEGM; Aschwanden et al., 2012)
% 2: modified enthalpy gradient model (MEGM; Hewitt I. and Schoof C., 2017)
% 3: isothermal

isFlowband = 1;
% true  (1): flowband mode
% false (0): flowline mode