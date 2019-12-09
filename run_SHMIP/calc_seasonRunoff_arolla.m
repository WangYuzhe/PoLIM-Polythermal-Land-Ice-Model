function seasonRunoff = calc_seasonRunoff_arolla(arrayTime, numTimeStep, deltaT)
% Date: 2019-5-17
% SHMIP Suite D

global hS M SPD SPY

DDF = 0.01/SPD; % degree dat factor [m K-1 s-1]
lr = -0.0065;   % lapse rate [K m-1]

Tref = -10*cos(2*pi*arrayTime/SPY) + deltaT; % Temperature at 0 m [degC]
delta_hS = hS - hS(end);
Tz = ones(numTimeStep,1)*delta_hS*lr + Tref*ones(1,M); % <numTimeStep * M>

seasonRunoff = DDF*Tz; % [m s-1]
seasonRunoff(seasonRunoff<0) = 0;

end    