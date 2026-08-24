% Date: 2019-5-17
% Date: 2019-11-10

clc
clearvars
clearvars -global

global SPD SPY M xi rhow g hB dx H

set_ice_parameters;

% ICE SHEET TOPOGRAPHY
ice_L = 6e3;
dx = 60;
xi = 0:dx:ice_L;
M = length(xi);

hS = 100*(xi+200).^(1/4) + xi/60 - (2*10^10)^(1/4) + 1;
hS = fliplr(hS);

% TIME SETTING
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

% WATER SOURCES
basalMelt = 2*5.79e-7*ones(1,M); % [m s-1]
mouInput = zeros(1,M);
seasonRunoff = zeros(numTimeStep,M);

% PARAMETERS
iter_threshold = 1;
ev = 1e-3; % englacial void fraction
ub = 1e-6*ones(1,M);
gammab = [0.05, 0, -0.1, -0.5, -0.7];
visc_ice_b = 1e12*ones(1,M);

for i = 1:length(gammab)
    fprintf('E%d is ruuning!\n', i)
    hB = (hS(end) - 6000*gammab(i)) / (6000^2) * xi.^2 + gammab(i) * xi;
    
    hB = fliplr(hB);    
    H = hS - hB + 0.01;
    
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
            
            if iTimeStep == 1
                phi_nm1 = rhow*g*hB;
                Hw_nm1 = 1e-1*ones(1,M); % [m]
                dHwdx_nm1 = [(Hw_nm1(2)-Hw_nm1(1))/dx,...
                    (Hw_nm1(3:end)-Hw_nm1(1:end-2))/(2*dx),...
                    (Hw_nm1(end)-Hw_nm1(end-1))/dx];
                dHw3dx_nm1 = 3*Hw_nm1.^2.*dHwdx_nm1;
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
    
    resultFileName = strcat('N_E', num2str(i));
    eval([resultFileName ' = Neff;'])
    eval(['save ' resultFileName ' ' resultFileName])    
end