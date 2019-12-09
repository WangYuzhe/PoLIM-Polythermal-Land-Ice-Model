clear, clc

% PARAMETERS
SPD = 24*3600; % [s d-1]
SPY = 365*SPD; % [s yr-1]

% GLACIER TOPOGRAPHY
ice_L = 6e3;
ice_W = 1e3;
dx = 60;
xi = 0:dx:ice_L;
M = length(xi);
hS = 100*(xi+200).^(1/4) + xi/60 - (2*10^10)^(1/4) + 1;
gammab = 0.05;
hB = (hS(end)-6000*gammab)/(6000^2)*xi.^2 + gammab*xi;

hS = fliplr(hS);
hB = fliplr(hB);

H = hS - hB + 0.01;

% TIME SETTING
dt_hydro = 1/24*SPD;
endTime = 1*SPY;
timeSpan = 0:dt_hydro:(endTime-dt_hydro);
numTimeStep = length(timeSpan);

% WATER SOURCES
basalMelt = 7.93e-11*ones(1,M); % [m s-1]
mouInput = zeros(1,M);

% PARAMETERS
iterThreshold = 4;
deltaT = [-6, -3, 0, 3, 6];

for i = 1:length(deltaT)
    fprintf('F%d is ruuning!\n', i)
    seasonRunoff = calc_seasonRunoff(hS, M, timeSpan, numTimeStep, deltaT(i));
    totalRechg = sum(seasonRunoff, 2)*dx*ice_W;
    
    resultFileName = strcat('totalRechg', '_F', num2str(i), ' ');
    eval([ 'save ' resultFileName ' totalRechg '])    

end