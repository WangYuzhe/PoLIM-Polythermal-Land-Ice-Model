clc
clearvars
clearvars -global

tic
global N M Ms hS hB H iTimeStep enth_lst

addpath('.\geo_inputs')

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params();
SPY = p.SPY;
iter_max = p.iter_max;
p.Hmin = 0.1;
p.layers = 31;
p.Qgeo = 0.04;

p.is_flowband = 1;
p.is_thk_evolv = 0;
p.is_surf_relax = 0;
p.type_BBC = 'CFlaw_polyT'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'; 'CFlaw_polyT'
p.type_enth_solver = 'MEGM';
p.is_auto_enth_BBC = 1;
p.type_Arrhenius = 'cuffey';
p.is_duval = 1;
p.type_valley = 'trapz';

%% GLACIER GEOMETRY
%
set_ice_geometry('./geo_inputs/geo_arolla.mat', p);

%% TIME SETTING
%
%
dt = 1; % time step for velocity solver [a]
startTime = 1; % start time
endTime = 10; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt, startTime, endTime);

%% INITIALIZATION
%
%
% Related to Ney-Glen law
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

% Related to velocity solver
At_u = zeros(N,M,numTimeStep);

% Related to enthalpy solver
% At_E = zeros(N,M,numTimeStep);
% At_T = zeros(N,M,numTimeStep);
% At_CTS = zeros(numTimeStep,M);
% At_m_basal = zeros(numTimeStep,M); % basal melt rate

% Related to glacier evolution
At_hS = zeros(numTimeStep,M);
At_hB = zeros(numTimeStep,M);
At_H = zeros(numTimeStep,M);
At_SMB = zeros(numTimeStep,M);

% SBC and IC for the enthalpy solver
if ~strcmp(p.type_enth_solver, 'isoT')
    Esbc = set_thermalSBC(p);
    Eini = get_initial_enthalpy(Esbc, p); % result_init.E;
end

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('Time step: %d, Year: %.1f \n', iTimeStep, arrayTime(iTimeStep))

    set_staggered_grid();

    % calculate the horizontal velocity
    [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_smedt(visc_s, visc, AGlen_s, p);
    % [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_pimentel(visc_s, visc, AGlen_s, p);
    
    % calculate the vertical velocity
    [w_vs, w] = get_ice_w(u_s, u);

    % calculate the enthalpy
    if strcmpi(p.type_enth_solver,'MEGM')
        % calculate friction heat
        frictionHeat = calc_frictionHeat(u_s, u, visc_s, p);

        % calculate ice enthalpy
        [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP, qw_TEMP_darcy] =...
            solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt, Esbc, Eini, p);

        % update the flow rate factor
        AGlen_s = get_AGlen(T, omega, CTS, p);

        % program termination
        if any(isnan(u(:))) || any(isnan(E(:)))
            break
        end
    end

    % print results
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))

    %% ICE THICKNESS EVOLUTION
    %
    %
    if p.is_thk_evolv
        Hn = H;
        %         SMB = 5e-3*(hS - 3150); % 3940 is for Arolla; 4000~4100 for incline
        hSn = hS;
        zref = 2500;
        SMB = calc_SMB(2500, -5.5, hSn);
        %         SMB = calc_SMB_test();

        SMB(Hn<p.Hmin) = 0;
        At_SMB(iTimeStep, :) = SMB;

        Hnp1 = get_evolution_continuity(Hn, u, SMB, dt);
        Hnp1(Hnp1<p.Hmin) = p.Hmin;
        hS = hB + Hnp1;
        H = Hnp1;
    end

    if p.is_surf_relax
        Hn = H;
        hSn = hS;
        SMB = zeros(1,M);
        SMB(Hn<p.Hmin) = 0;
        At_SMB(iTimeStep, :) = SMB;

        hSnp1 = get_evolution_kinematic(hSn, u, w, SMB, dt);
        Hnp1 = hSnp1-hB;
        Hnp1(Hnp1<p.Hmin) = p.Hmin;
        hS = hB + Hnp1;
        H = Hnp1;
    end
    %%
    %
    % STORAGE RESULTS
    At_u(:,:,iTimeStep) = u;
    At_hS(iTimeStep, :) = hS;
    At_hB(iTimeStep, :) = hB;
    At_H(iTimeStep, :) = H;

    if ~strcmpi(p.type_enth_solver,'isoT')
        % STORAGE RESULTS
        enth_lst.E = E;
        enth_lst.omega = omega;
        enth_lst.Kappa_vs = Kappa_vs;
        enth_lst.thk_w = thk_w;
        enth_lst.thk_TEMP = thk_TEMP;
        enth_lst.is_TEMP = is_TEMP;
        enth_lst.qw_TEMP_darcy = qw_TEMP_darcy;
        enth_lst.CTS = CTS;
        %
        %         At_E(:,:,iTimeStep) = E;
        %         At_T(:,:,iTimeStep) = T;
        %         At_omega(:,:,iTimeStep) = omega;
        %         At_CTS(iTimeStep,:) = CTS;
        %         At_m_basal(iTimeStep,:) = m_basal;
        %
    end
end
toc

% plots

figure(1)
plot_uField_uSurf
if  ~strcmpi(p.type_enth_solver,'isoT')
    figure(2)
    plot_enthalpy_field
end

% plot_along_flow_temp_profiles

% result_init.AGlen_s = AGlen_s;
% result_init.E = E;
% save('result_warm.mat', 'result_init');

