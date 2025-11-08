clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms hS hB H iTimeStep

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

p.type_BBC = 'zero'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'
p.lambda_max = 2; % wavelength of the dominant bedrock bumps [m]
p.m_max = 0.4; % maximum slope of the dominant bedrock bumps []; ref: 0.5

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

visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

% Related to velocity solver
At_u = zeros(N,M,numTimeStep);

% Related to glacier evolution
At_hS = zeros(numTimeStep,M);
At_hB = zeros(numTimeStep,M);
At_H = zeros(numTimeStep,M);
At_SMB = zeros(numTimeStep,M);

%% MAIN
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d, year: %d \n', iTimeStep, arrayTime(iTimeStep))

    set_staggered_grid();

    [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_smedt(visc_s, visc, AGlen_s, p);
    % [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_pimentel(visc_s, visc, AGlen_s, p);
    %% ICE THICKNESS EVOLUTION
    %
    %
    if p.is_thk_evolv
        Hn = H;
        hSn = hS;

        zref = hSn(end);
        m_b = calc_SMB(zref, -2, hSn);

        m_b(Hn<p.Hmin) = 0;
        At_SMB(iTimeStep, :) = m_b;

        Hnp1 = get_evolution_continuity(Hn, u, m_b, dt_u);
        Hnp1(Hnp1<=p.Hmin) = p.Hmin;
        hS = hB + Hnp1;
        H = Hnp1;
    end

    if p.is_surf_relax
        Hn = H;
        hSn = hS;
        SMB = zeros(1,M);
        SMB(Hn<p.Hmin) = 0;
        At_SMB(iTimeStep, :) = SMB;

        hSnp1 = get_evolution_kinematic(hSn, u, w, SMB, dt_u);
        Hnp1 = hSnp1-hB;
        Hnp1(Hnp1<p.Hmin) = p.Hmin;
        hS = hB + Hnp1;
        H = Hnp1;
    end
    %%
    %
    % RESULTS
    At_u(:,:,iTimeStep) = u;
    At_hS(iTimeStep, :) = hS;
    At_hB(iTimeStep, :) = hB;
    At_H(iTimeStep, :) = H;

    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))

    %%
    if sum(isnan(u(:))) || sum(H<=p.Hmin)==M
        break
    end

end
toc

% plot the results
% plot_uField_uSurf

