function set_ice_parameters_ismip_hom(id_type_BBC, id_type_geometry,...
    id_type_valley, id_isFlowband)
global SPY rho g n de0 epsilon iter_max lambda_max m_max epr...
    type_BBC type_valley isFlowband type_geometry 

% fixed physical constants
SPY = 31556926;    % year is this many seconds (i.e. 365.2422 days)
rho = 910;         % density of ice [kg m-3]
g = 9.81;          % accel of gravity [m s-2]
de0 = 1e-30;       % small number in case of singularity [yr-2]
n = 3;             % exponent in Glen's flow law []

% adjustable parameters
epsilon = 0;       % enhanced stress-free condition
iter_max = 50;     % iterations

lambda_max = 4; % ref: 4
m_max = 0.3; % ref: 0.3
epr = 1; % N=Pi-Pw, epr=0, Pw=Pi; epr=1, Pw=0.

type_BBC = id_type_BBC;
% 1: no-slip bed
% 2: Coulomb friction law
% 3: linear friction law
% 4: linear friction law for ISMIP-HOM E2
% 5: simplified linear friction law

type_geometry = id_type_geometry;
% 1: Haut d'Arolla
% 2: Haut d'Arolla (resampled to 250 grid points)
% 3: Kleiner Exp. A (parallel-sided slab)
% 4: Kleiner Exp. B (polythermal parallel-sided slab)
% 5: Hewitt-Schoof Ice Cap Experiment

type_valley = id_type_valley;
% 1: trapezoid
% 2: Svessen: y=ax^b
% 3: rectangular

isFlowband = id_isFlowband;
% true  (1): flowband mode
% false (0): flowline mode
end