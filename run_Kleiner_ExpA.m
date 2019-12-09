% Date: 2017-7-23
clc
clearvars
clearvars -global

tic
global N M Cp Tref At_E At_T At_isTemperate iTimeStep At_Hw At_Ht At_Kappa_s At_omega

set_ice_parameters_Kleiner_ExpA;
set_ice_geometry();

% Time setting
dt = 100; % [a]
[arrayTime, numTimeStep] = set_time_step(dt, 300001);

% Initialization
%    Related to enthalpy solver
At_E                       = zeros(N, M, numTimeStep);
At_T                       = zeros(N, M, numTimeStep);
At_omega                   = zeros(N, M, numTimeStep);
At_Kappa_s                 = zeros(N-1, M, numTimeStep);
At_CTS                     = zeros(numTimeStep, M);
At_Ht                      = zeros(numTimeStep, M);
At_basalMeltRate           = zeros(numTimeStep, M);
At_Hw                      = zeros(numTimeStep, M);
At_basalT                  = zeros(numTimeStep, M);
At_isTemperate             = zeros(numTimeStep, M);
At_temperateWaterFlux      = zeros(N-1, M, numTimeStep);
At_drainToBed              = zeros(numTimeStep, M);

%    Velocity field
u = zeros(N,M);
u_s = zeros(N,M+1);
w = zeros(N,M);
w_vs = zeros(N,M);
strainHeat = zeros(N,M);

% Initial enthalpy field
Eini = Cp*(-30 + 273.15 - Tref)*ones(N,M);

% Time-dependent thermal surface boundary condition
Tsbc = [-30*ones(1,100000), -5*ones(1,50000), -30*ones(1,150001)] + 273.15;
Esbc = Cp*(Tsbc - Tref);

% options for the thermal model
is_auto_thermalBasalBC = 1;
type_thermalBasalBC = 1;
has_Greve_drainage = 1;

for iTimeStep = 1:numTimeStep
    yr = arrayTime(iTimeStep);
    
    % solve the enthalpy balance equation
    [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, drainToBed] = ...
        solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc(yr), Eini, is_auto_thermalBasalBC, type_thermalBasalBC, has_Greve_drainage);

    At_E(:,:,iTimeStep)                  = E;
    At_T(:,:,iTimeStep)                  = T;
    At_omega(:,:,iTimeStep)              = omega;    
    At_Kappa_s(:,:,iTimeStep)            = Kappa_s;
    At_CTS(iTimeStep,:)                  = CTS;
    At_Ht(iTimeStep,:)                   = Ht;
    At_basalMeltRate(iTimeStep,:)        = basalMeltRate;
    At_basalT(iTimeStep,:)               = T(1,:);
    At_temperateWaterFlux(:,:,iTimeStep) = temperateWaterFlux;
    At_drainToBed(iTimeStep,:)           = drainToBed;
    
    % calculate the basal water layer thickness
    if iTimeStep == 1
        At_Hw(iTimeStep,:) = zeros(1, M) + dt*basalMeltRate;
    else
        At_Hw(iTimeStep,:) = At_Hw(iTimeStep-1,:) + dt*basalMeltRate;
    end
    
    % logic value for basal thermal state
    logic1 = (At_Hw(iTimeStep,:)>0 & At_isTemperate(iTimeStep,:)==0);
    At_isTemperate(iTimeStep,:) = At_isTemperate(iTimeStep,:) | logic1;
end
toc

plot_enthalpy_ExpA

save result_enthA arrayTime At_basalT At_basalMeltRate At_Hw