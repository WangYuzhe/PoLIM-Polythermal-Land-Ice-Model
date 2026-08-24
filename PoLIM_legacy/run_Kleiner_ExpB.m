% Date: 2017-7-23
clc
clearvars
clearvars -global

tic
global N M zeta H rho g Cp Tref AGlen iTimeStep At_E At_T At_omega At_Kappa_s...
    At_Hw At_Ht SPY At_isTemperate

set_ice_parameters_Kleiner_ExpB;
set_ice_geometry();

% Time setting
% CFL condition: |vz*dt/dz|<=1 ---> dt<=|dz/vz|
% N=21, dz=10 m, dt <= 50; N=401, dz=0.5 m, dt<=2.5
dt = 1; % [a];
[arrayTime, numTimeStep] = set_time_step(dt, 2000);

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

%    Velocity field
alpha = 4*pi/180; % inclination angle
vz = -0.2; % [m a-1]

u = zeros(N,M);
u_s = zeros(N,M+1);
w = vz*ones(N,M);
w_vs = vz*ones(N,M);
strainHeat = 2*AGlen*(rho*g*sin(alpha))^4*H^4*(1 - zeta).^4*SPY;

% Initial enthalpy field
Eini = Cp*(-1.5+273.15- Tref)*ones(N,M);

% Thermal surface boundary condition
Tsbc = (-3 + 273.15)*ones(N,M);
Esbc = Cp*(Tsbc - Tref);

% options for the thermal model
is_auto_thermalBasalBC = 1;
type_thermalBasalBC = 1;
has_Greve_drainage = 0;

for iTimeStep = 1:numTimeStep   
    % solve the enthalpy balance equation
    [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux] =...
        solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc, Eini, is_auto_thermalBasalBC, type_thermalBasalBC, has_Greve_drainage);
    
    At_E(:,:,iTimeStep)                  = E;
    At_T(:,:,iTimeStep)                  = T;
    At_omega(:,:,iTimeStep)              = omega;    
    At_Kappa_s(:,:,iTimeStep)            = Kappa_s;
    At_CTS(iTimeStep,:)                  = CTS;
    At_Ht(iTimeStep,:)                   = Ht;
    At_basalMeltRate(iTimeStep,:)        = basalMeltRate;
    At_basalT(iTimeStep,:)               = T(1,:);
    At_temperateWaterFlux(:,:,iTimeStep) = temperateWaterFlux;
    
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

result = [zeta, At_E(:,:,end), At_T(:,:,end), At_omega(:,:,end)];
plot_enthalpy_ExpB

if N == 21
    save result_enthB_10m result
else
    save result_enthB_0p5m result
end