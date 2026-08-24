% Run the ice cap experiment in Hewitt&Schoof (2017)
% schemeID = 1: the enthalpy solver uses the modified enthalpy gradient
% method

% schemeID = 2: the enthalpy solver uses the standard enthalpy gradient
% method without drainage

% schemeID = 3: the enthalpy solver uses the standard enthalpy gradient
% method with drainage

% Date: 2018-4-15
% Author: Wang Yuzhe
% Enthalpy is defined as E = Cp*(T - Tref)

clc
clearvars
clearvars -global

tic
global M N zeta hB H AGlen rho g Cp Tref SPY iTimeStep dhSdx
global At_E At_T At_omega At_Kappa_s At_Hw At_Ht At_temperateWaterFluxDarcy

set_ice_parameters_Hewitt;
set_ice_geometry();

% Time setting
dt = 1; % [a]
[arrayTime, numTimeStep] = set_time_step(dt, 4000);

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
At_temperateWaterFluxDarcy = zeros(N-1,M,numTimeStep);
At_drainToBed              = zeros(numTimeStep, M);

%-----------------solve the velocity field-----------------
set_staggered_grid();

u = zeros(N, M);
strainHeat = zeros(N, M);
for i = 1:M
    for j = 1:N
        u(j,i) = -2*AGlen*(rho*g*dhSdx(i)).^3/4*(H(i).^4 - ((1-zeta(j))*H(i)-hB(i)).^4);
        strainHeat(j,i) = 2*AGlen*(rho*g*dhSdx(i)).^4*((1-zeta(j)).*H(i)-hB(i)).^4;
    end
end
u_s = main2staggerX(u);
[w_vs, w] = get_ice_w(u_s, u);

% convert [m s-1] to [m a-1]
u = u*SPY;
u_s = u_s*SPY;
w = w*SPY;
w_vs = w_vs*SPY;
strainHeat = strainHeat*SPY;
%-----------------solve the velocity field-----------------

% Thermal surface boundary condition
Ts = -1;
Esbc = Cp*(Ts + 273.15- Tref)*ones(1, M);

% Initial enthalpy field
Eini = Cp*(262.15 - Tref)*ones(N, M);

schemeID = 3;
% 1: MEGM
% 2: SEGM0 (without drainage)
% 3: SEGM1 (with drainage)

for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    % solve the enthalpy balance equation
    if schemeID == 1
        [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, temperateWaterFluxDarcy] =...
            solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc, Eini, 0, 2);

    elseif schemeID == 2
        [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, drainToBed] = ...
            solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc, Eini, 0, 2, 0);
    else
        [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, drainToBed] = ...
            solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc, Eini, 0, 2, 1);
    end
    
    At_E(:,:,iTimeStep)                  = E;
    At_T(:,:,iTimeStep)                  = T;
    At_omega(:,:,iTimeStep)              = omega;
    At_Kappa_s(:,:,iTimeStep)            = Kappa_s;
    At_CTS(iTimeStep,:)                  = CTS;
    At_Ht(iTimeStep,:)                   = Ht;
    At_basalMeltRate(iTimeStep,:)        = basalMeltRate;
    At_basalT(iTimeStep,:)               = T(1,:);
    At_temperateWaterFlux(:,:,iTimeStep) = temperateWaterFlux;
    
    if schemeID==1        
        At_temperateWaterFluxDarcy(:,:,iTimeStep) = temperateWaterFluxDarcy;
    end

    if schemeID ~= 1
        At_drainToBed(iTimeStep,:) = drainToBed;
    end
    
    % calculate the basal water layer thickness
    if iTimeStep == 1
        At_Hw(iTimeStep,:) = zeros(1,M) + dt*basalMeltRate;
    else
        At_Hw(iTimeStep,:) = At_Hw(iTimeStep-1,:) + dt*basalMeltRate;
    end
    
    % logic value for basal thermal state
    logic1 = (At_Hw(iTimeStep,:)>0 & At_isTemperate(iTimeStep,:)==0);
    At_isTemperate(iTimeStep,:) = At_isTemperate(iTimeStep,:) | logic1;
end
toc

if schemeID == 1
    plot_field_MEGM
    if Ts == -10
        save result_icecap_MEGM_m10 E T omega temperateWaterFlux CTS xx yy xi hS hB u
    else
        save result_icecap_MEGM_m1 E T omega temperateWaterFlux CTS xx yy xi hS hB u
    end
elseif schemeID == 2
    plot_field_SEGM0
    if Ts == -10
        save result_icecap_SEGM0_m10 E T omega temperateWaterFlux CTS xx yy xi hS hB u
    else
        save result_icecap_SEGM0_m1 E T omega temperateWaterFlux CTS xx yy xi hS hB u
    end
else
    plot_field_SEGM1
    if Ts == -10
        save result_icecap_SEGM1_m10 E T omega temperateWaterFlux CTS xx yy xi hS hB u drainToBed
    else
        save result_icecap_SEGM1_m1 E T omega temperateWaterFlux CTS xx yy xi hS hB u drainToBed
    end    
end