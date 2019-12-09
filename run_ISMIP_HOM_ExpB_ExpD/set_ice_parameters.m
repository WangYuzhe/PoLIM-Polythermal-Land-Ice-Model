global SPY rho g n de0 AGlen iter_max experiment

SPY = 31556926;    % year is this many seconds (i.e. 365.2422 days)
rho = 910;         % density of ice [kg m-3]
g = 9.81;          % accel of gravity [m s-2]
n = 3;             % exponent in Glen's flow law []
AGlen = 1e-16;     % [Pa-3 a-1]

de0 = 1e-30;       % small number in case of singularity [yr-2]
iter_max = 600;     % iterations

experiment = 2;
% 1: Experiment B
% 2: Experiment D