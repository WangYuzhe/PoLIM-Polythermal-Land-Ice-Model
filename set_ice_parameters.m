global SPY rho g n de0 AGlen
global epsilon iter_max type_BBC type_valley lambda_max m_max epr isFlowband

% fixed physical constants
SPY = 31556926;    % year is this many seconds (i.e. 365.2422 days)
rho = 910;         % density of ice [kg m-3]
g = 9.81;          % accel of gravity [m s-2]
de0 = 1e-30;       % small number in case of singularity [yr-2]
n = 3;             % exponent in Glen's flow law []
AGlen = 1e-16; % [Pa-3 a-1]

% adjustable parameters
epsilon = 0;       % enhanced stress-free condition
iter_max = 50;     % iterations

type_BBC = 3;
% 1: no-slip bed
% 2: Coulomb friction law
% 3: linear friction law

type_valley = 3;
% 1: trapezoid
% 2: Svessen: y=ax^b
% 3: rectangular

isFlowband = 0;
% true  (1): flowband mode
% false (0): flowline mode

lambda_max = 4; % ref: 4
m_max = 0.3; % ref: 0.3
epr = 1; % N=Pi-Pw, epr=0, Pw=Pi; epr=1, Pw=0.