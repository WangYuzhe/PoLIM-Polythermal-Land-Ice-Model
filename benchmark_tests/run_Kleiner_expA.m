clc
clearvars
clearvars -global

addpath('..', '-frozen')

tic
global N M dzetadx iTimeStep enth_lst

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
% p = params_hewitt_schoof_slab();
p = params_Kleiner_ExpA();
SPY = p.SPY;

p.layers = 101;
%% GLACIER GEOMETRY
%
%
set_ice_geometry('../geo_inputs/geo_Kleiner_ExpA.mat', p);
dzetadx = zeros(N,M);
%% TIME SETTING
%
%
dt = 100; % time step for velocity solver [a]
[arrayTime, numTimeStep] = set_time_step(dt,1,300001);

%% INITIALIZATION
%
%
u = zeros(N,M);
u_s = zeros(N,M+1);
w = zeros(N,M);
w_vs = zeros(N,M);
strainHeat = zeros(N,M);
frictionHeat = zeros(1,M);

%  Related to enthalpy solver
At_T = zeros(N,M,numTimeStep);
At_is_TEMP = zeros(numTimeStep,M);
At_thk_w = zeros(numTimeStep,M); % water thickness
At_m_basal = zeros(numTimeStep,M); % basal melt rate

% SBC and IC for the enthalpy solver
Tsbc = [-30*ones(1,100000), -5*ones(1,50000), -30*ones(1,150001)] + 273.15;
Eini = p.Cp*(-30 + 273.15 - p.Tref)*ones(N,M);

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    Esbc = p.Cp*(Tsbc(trueTime)-p.Tref)*ones(N,M);
    
    % solver SEGM
    [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP_drain_greve, qw_TEMP_diffu] = ...
        solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt, Esbc, Eini, p);

    % STORAGE RESULTS
    enth_lst.E = E;
    enth_lst.omega = omega;
    enth_lst.Kappa_vs = Kappa_vs;
    enth_lst.thk_w = thk_w;
    enth_lst.thk_TEMP = thk_TEMP;
    enth_lst.is_TEMP = is_TEMP;
        
    At_T(:,:,iTimeStep) = T;
    At_m_basal(iTimeStep,:) = m_basal;
    At_thk_w(iTimeStep,:) = thk_w;
    At_is_TEMP(iTimeStep,:) = is_TEMP;
    
end
toc

% figure
load('enthA_analy_result.mat')
plot_benchmark.plot_enthalpy_ExpA(arrayTime,At_T,At_m_basal,At_thk_w,At_is_TEMP,enthA_analy_result)