% Date: 2019-6-3

clear,clc

% PARAMETERS
SPD = 24*3600; % [s d-1]
SPY = 365*SPD; % [s yr-1]

% ICE SHEET TOPOGRAPHY
ice_L = 100e3;
ice_W = 20e3;
dx = 1000;
xi = 0:dx:ice_L;
M = length(xi);
hS = 6*((xi+5000).^(1/2)-sqrt(5000))+1;
hB = zeros(1,M);

hS = fliplr(hS);
hB = fliplr(hB);

H = hS - hB + 0.01;

% index for the three bands
Hindex1 = find(xi==10000); Hindex2 = find(xi==15000); % highest band
Mindex1 = find(xi==50000); Mindex2 = find(xi==55000); % middle band
Lindex1 = find(xi==85000); Lindex2 = find(xi==90000); % lower band

% TIME SETTING
dt_hydro = 1/24*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
timeSpan = 0:dt_hydro:(endTime-dt_hydro);
numTimeStep = length(timeSpan);

% PARAMETERS
deltaT = [-4, -2, 0, 2, 4];

for i = 1:length(deltaT) %
    fprintf('D%d is running!\n', i)
    seasonRunoff = calc_seasonRunoff(hS, M, timeSpan, numTimeStep, deltaT(i));
    totalRechg = sum(seasonRunoff, 2)*ice_W*ice_L;
       
    bandRechg = [sum(seasonRunoff(:, 2:Hindex1),2)*dx*ice_W, sum(seasonRunoff(:, 2:Mindex1),2)*dx*ice_W,...
        sum(seasonRunoff(:, 2:Lindex1),2)*dx*ice_W]; % upper, middle, lower
    
    resultFileName = strcat('Rechg', '_D', num2str(i), ' ');
    eval([ 'save ' resultFileName ' totalRechg ' 'bandRechg'])    
end