% Date: 2019-5-17
% Date: 2019-11-10

clc
clearvars
clearvars -global

global SPD SPY type_glacier_geometry M xi rhow g hB dx

set_ice_parameters;

% ice sheet case is not good
type_glacier_geometry = 2;
% 1: haut d'arolla
% 2: ice sheet
% 3: valley glacier

set_ice_geometry();

% TIME SETTING
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

% WATER SOURCES
basalMelt_range = [7.93e-11, 1.59e-9, 5.79e-9, 2.5e-8, 4.5e-8, 5.79e-7]; % [m s-1]
mouInput = zeros(1,M);
seasonRunoff = zeros(numTimeStep,M);

% PARAMETERS
iter_threshold = 1;
ev = 0; % englacial void fraction
ub = 1e-6*ones(1,M);
visc_ice_b = 1e12*ones(1,M);

for i = 1:length(basalMelt_range)
    fprintf('A%d is running!\n', i)
    basalMelt = basalMelt_range(i)*ones(1,M);
    
    % INITIAL CONDITIONS
    At_Hw = zeros(numTimeStep, M);
    At_dHw3dx = zeros(numTimeStep,M);
    At_phi = zeros(numTimeStep,M);
    At_N = zeros(numTimeStep, M);
    iter = 0;    
    while 1
        iter = iter + 1;
        for iTimeStep = 1:numTimeStep
            fprintf('Time Step: %d is running!\n', iTimeStep)
            
            if iTimeStep == 1 && iter ==1
                phi_nm1 = rhow*g*hB;                
                Hw_nm1 = 1e-1*ones(1,M); % [m]
                dHwdx_nm1 = [(Hw_nm1(2)-Hw_nm1(1))/dx,...
                    (Hw_nm1(3:end)-Hw_nm1(1:end-2))/(2*dx),...
                    (Hw_nm1(end)-Hw_nm1(end-1))/dx];
                dHw3dx_nm1 = 3*Hw_nm1.^2.*dHwdx_nm1;
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
        end
        if iter > iter_threshold
            break
        end
    end
    
    resultFileName = strcat('N_A', num2str(i));
    eval([resultFileName ' = Neff;'])
    eval(['save ' resultFileName ' ' resultFileName])    
end