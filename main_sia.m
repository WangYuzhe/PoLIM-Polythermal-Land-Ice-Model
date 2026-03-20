clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms u_s_lst iTimeStep

%% PHYSICAL CONSTANTS AND PARAMETERS
%
p = params();
SPY = p.SPY;

p.Hmin = 0.1;
p.layers = 31;
p.type_valley = 'rect';

p.is_flowband = 0;
p.is_thk_evolv = 0;
p.is_surf_relax = 0;

p.type_BBC = 'CFlaw_isoT'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'
p.lambda_max = 4; % wavelength of the dominant bedrock bumps [m]
p.m_max = 0.3; % maximum slope of the dominant bedrock bumps []

%% GLACIER GEOMETRY
%
set_ice_geometry('./geo_inputs/geo_arolla.mat', p);

%% TIME SETTING
%
dt_u = 1.0; % time step for velocity solver [a]
startTime = 1; % start time of model run [a]
endTime = 1; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, startTime, endTime);

%% INITIALIZATION
%
% Related to Ney-Glen law
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]
visc_s_init = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc_init = zeros(N,M) + 1e13/SPY; % [Pa a]

% data_init = load('./outinit_large_slide.mat');
% visc_s_init = data_init.visc_s; % [Pa a]
% visc_init = data_init.visc; % [Pa a]
% u_s_lst = data_init.u_s;

% outinit = './outinit_large_slide.mat';

%% SOLVER
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d, year: %d \n', iTimeStep, arrayTime(iTimeStep))

    set_staggered_grid();
    % if iTimeStep==1
    %     u_s_lst = sia_u(AGlen_s, p);
    % end
    
    % [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_smedt_test(visc_s_init,visc_init,AGlen_s,p);
    [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_smedt(visc_s_init,visc_init,AGlen_s,p);

    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
end
toc

%% PLOT
%
figure()
plot_uField_uSurf

%% OUTPUT
%
% save(outinit, 'visc_s', 'visc', 'u_s', 'u')
