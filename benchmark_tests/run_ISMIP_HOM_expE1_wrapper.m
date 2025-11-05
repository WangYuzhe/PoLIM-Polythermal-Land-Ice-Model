clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms iTimeStep

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params();
SPY = p.SPY;
iter_max = p.iter_max;
p.Hmin = 0.01;
p.layers = 31;
p.type_valley = 'rect';

p.is_flowband = 0;
p.type_BBC = 'zero'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'; 'CFlaw_polyT'

p.type_Arrhenius = 'cuffey';
p.type_enth_solver = 'isoT';

%% GLACIER GEOMETRY
%
set_ice_geometry('../geo_inputs/geo_arolla', p);

%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 1; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, 1, endTime);

%% INITIALIZATION
%
%
% Related to Ney-Glen law
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]

visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)

    trueTime = arrayTime(iTimeStep);

    set_staggered_grid();
    [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_chatgpt(visc_s, visc, AGlen_s, p);

    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))

end
toc
plot_uField_uSurf


global xi
normalized_xi = (xi - xi(1))./(xi(end) - xi(1));
[tauxz, ~] = calc_tauxz(u_s, visc_s);

figure()
set(gcf, 'position', [146,332,1000,350]);
subplot(1,2,1)
load('ismip_hom_us_mean_std.mat')
plot_benchmark.plot_shadedErrorBar(NFS_e000, FS_e000, 'Velocity (m a^{-1})')
hold on
plot(normalized_xi, u(end,:), 'k--', 'linewidth', 1)

subplot(1,2,2)
load('ismip_hom_expE_taub_mean_std.mat')
plot_benchmark.plot_shadedErrorBar(NFS_e000, FS_e000, 'Basal shear stress (kPa)')
hold on
plot(normalized_xi, tauxz/1e3, 'k--', 'linewidth', 1)