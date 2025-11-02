clc
clearvars
clearvars -global

tic

addpath('..')

global N M H dzetadx zeta iTimeStep enth_lst

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params_Kleiner_ExpB();
SPY = p.SPY;

p.layers = 201;
%% GLACIER GEOMETRY
%
%
set_ice_geometry('../geo_inputs/geo_kleiner_expB', p);
dzetadx = zeros(N,M);
%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
start_time = 1;
end_time = 2000; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, start_time, end_time);

%% INITIALIZATION
%
%
alpha = 4*pi/180; % inclination angle
vz = -0.2; % [m a-1] vertical velocity

u = zeros(N,M);
u_s = zeros(N,M+1);
w = vz*ones(N,M);
w_vs = vz*ones(N,M);
strainHeat = 2*p.AGlen*(p.rhoi*p.g*sin(alpha)).^4.*H.^4.*(1-zeta).^4*SPY;
frictionHeat = zeros(1,M);

%  Related to enthalpy solver
At_E = zeros(N,M,numTimeStep);
At_T = zeros(N,M,numTimeStep);
At_omega = zeros(N,M,numTimeStep);
At_Kappa_vs = zeros(N-1,M,numTimeStep);
At_CTS = zeros(numTimeStep,M);
At_is_TEMP = zeros(numTimeStep,M);
At_thk_TEMP = zeros(numTimeStep,M);
At_thk_w = zeros(numTimeStep,M); % water thickness
At_m_basal = zeros(numTimeStep,M); % basal melt rate

% SBC and IC for the enthalpy solver
Esbc = p.Cp*(-3 + 273.15 - p.Tref);
Eini = p.Cp*(-1.5 + 273.15 - p.Tref)*ones(N,M);
%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    % calculate the enthalpy
    [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP_drain_Greve, qw_TEMP_diffu] = ...
        solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
    
    % STORAGE RESULTS
    enth_lst.E = E;
    enth_lst.omega = omega;
    enth_lst.Kappa_vs = Kappa_vs;
    enth_lst.thk_w = thk_w;
    enth_lst.thk_TEMP = thk_TEMP;
    enth_lst.is_TEMP = is_TEMP;
    
%     enth_lst.qw_TEMP_darcy = qw_TEMP_darcy;
    
end
toc

% figure
% plot_enthalpy_ExpB
load('enthB_analy_result.mat')
plot_benchmark.plot_enthalpy_ExpB(zeta,E,T,omega,enthB_analy_z,...
    enthB_analy_E,enthB_analy_T,enthB_analy_omega)