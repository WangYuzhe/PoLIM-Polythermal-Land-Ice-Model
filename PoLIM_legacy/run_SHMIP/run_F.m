% Date: 2019-5-17
% Date: 2019-11-10
clc
clearvars
clearvars -global

global SPD SPY type_glacier_geometry M xi rhow g hB dx

set_ice_parameters;

% ice sheet case is not good
type_glacier_geometry = 3;
% 1: haut d'arolla
% 2: ice sheet
% 3: valley glacier

ice_L = 6e3;
ice_W = 1e3;
set_ice_geometry();

% index for the three bands
Hindex1 = find(xi==600); Hindex2 = find(xi==900); % highest band
Mindex1 = find(xi==2700); Mindex2 = find(xi==3000); % middle band
Lindex1 = find(xi==5100); Lindex2 = find(xi==5400); % lower band

% TIME SETTING
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

% WATER SOURCES
basalMelt = 7.93e-11*ones(1,M); % [m s-1]
mouInput = zeros(1,M);
deltaT = [-6, -3, 0, 3, 6];

% PARAMETERS
iter_threshold = 50;
ev = 1e-3; % englacial void fraction
ub = 1e-6*ones(1,M);
visc_ice_b = 1e12*ones(1,M);

for i = 1:length(deltaT)
    fprintf('D%d is running!\n', i)
    seasonRunoff = calc_seasonRunoff(arrayTime, numTimeStep, deltaT(i));
    totalRechg = sum(seasonRunoff, 2)*ice_W*ice_L;
    
    % INITIAL CONDITIONS
    At_Hw = zeros(numTimeStep, M);
    At_dHw3dx = zeros(numTimeStep,M);
    At_phi = zeros(numTimeStep,M);
    At_N = zeros(numTimeStep, M);
    At_qw = zeros(numTimeStep, M);
    iter = 0;
    while 1
        iter = iter + 1;
        for iTimeStep = 1:numTimeStep
            fprintf('Time Step: %d is running!\n', iTimeStep)
            
            if iTimeStep == 1 && iter ==1
                load result_A1.mat
                phi_nm1 = phi_A1;
                Hw_nm1 = Hw_A1; % [m]
                dHw3dx_nm1 = dHw3dx_A1;
            elseif iTimeStep == 1 && iter > 1
                Hw_nm1 = At_Hw(end,:);
                dHw3dx_nm1 = At_dHw3dx(end,:);
                phi_nm1 = At_phi(end,:);
            else
                Hw_nm1 = At_Hw(iTimeStep-1,:);
                dHw3dx_nm1 = At_dHw3dx(iTimeStep-1,:);
                phi_nm1 = At_phi(iTimeStep-1,:);
            end
            
            [Neff, phi, qw, Hw, dHw3dx] = subHydro_storage_Hoffman(dt,...
                basalMelt, mouInput, seasonRunoff(iTimeStep,:), Hw_nm1, dHw3dx_nm1, phi_nm1, ev, ub, visc_ice_b);
            
            At_Hw(iTimeStep, :) = Hw;
            At_dHw3dx(iTimeStep, :) = dHw3dx;
            At_phi(iTimeStep, :) = phi;
            At_N(iTimeStep, :) = Neff;
            At_qw(iTimeStep, :) = qw;
        end
        if iter > iter_threshold
            break
        end
    end
    % Discharge for the band
    At_Qw = [mean(At_qw(:, Hindex1:Hindex2),2)*ice_W,...
        mean(At_qw(:, Mindex1:Mindex2),2)*ice_W,...
        mean(At_qw(:, Lindex1:Lindex2),2)*ice_W]; % upper, middle, lower
    
    % SAVE RESULTS
    resultFileName = strcat('result_F', num2str(i), ' ');
    eval([ 'save ' resultFileName ' At_N ' 'At_qw ' 'At_Qw ' 'totalRechg'])
end