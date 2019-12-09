function seasonRunoff = calc_seasonRunoff(arrayTime, numTimeStep, deltaT)
% Date: 2019-5-17
% SHMIP Suite D

global hS M SPD SPY

DDF = 0.01/SPD; % degree dat factor [m K-1 s-1]
lapseRate = -0.0075;   % lapse rate [K m-1]

T0m = -16.0*cos(2*pi*arrayTime/SPY) - 5 + deltaT; % Temperature at 0 m [degC]

seasonRunoff = DDF*(ones(numTimeStep,1)*hS*lapseRate + T0m*ones(1,M)); % [m s-1]
seasonRunoff(seasonRunoff<0) = 0;

end    