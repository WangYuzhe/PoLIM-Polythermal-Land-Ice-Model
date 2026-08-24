clc
clearvars
clearvars -global

tic
global N M Ms dx hS hB H iTimeStep enth_lst is_active

addpath('.\geo_inputs')
addpath('.\mini_tools')

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params();
SPY = p.SPY;
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
is_active = H > p.Hmin;

%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
dt_E = 0.1; % time step for enthalpy solver [a]

startTime = 1; % start time
endTime = 20; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, startTime, endTime);

p.CFL = 0.4;
p.dt_H_max = dt_u;
p.dt_H_min = 0.001;

%% SMB coupling
smb_driver = struct();
smb_driver.mode = 'auto';
smb_driver.fallback = 'fujita';
smb_driver.climate_file = 'gfdl_hist.mat';
smb_driver.calendar_start_year = 2008;
smb_driver.model_start_time = startTime;
smb_driver.hydro_start_month = 10;
smb_driver.zref = 2943;
smb_driver.fujita_zref = 2500;
smb_driver.Tma = -5.5;
smb_driver.state = struct();
smb_driver.last = struct();

n_sub = round(dt_u / dt_E);

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

    if ~any(is_active)
        fprintf('Glacier is extinct; stopping the simulation.\n')
        break
    end

    set_staggered_grid();

    % calculate the horizontal velocity
    [u_s, u, visc_s, visc, strainHeat] = solve_u_smedt(visc_s, visc, AGlen_s, p);
    % [u_s, u, visc_s, visc, strainHeat] = solve_u_pimentel(visc_s, visc, AGlen_s, p);

    % calculate the vertical velocity
    [w_vs, w] = get_ice_w(u_s, u);

    % calculate the enthalpy
    if strcmpi(p.type_enth_solver,'MEGM')
        % calculate friction heat
        frictionHeat = calc_frictionHeat(u_s, u, visc_s, p);

        for i_sub = 1:n_sub
            % calculate ice enthalpy
            [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP, qw_TEMP_darcy] =...
                solve_enth_MEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_E, Esbc, Eini, p);
        end
        % update the flow rate factor
        AGlen_s = get_AGlen(T, omega, CTS, p);

        % Publish the current thermal state before any geometry-driven
        % velocity recomputation during thickness subcycling.
        enth_lst.E = E;
        enth_lst.Kappa_vs = Kappa_vs;
        enth_lst.thk_w = thk_w;
        enth_lst.thk_TEMP = thk_TEMP;
        enth_lst.is_TEMP = is_TEMP;
        enth_lst.qw_TEMP_darcy = qw_TEMP_darcy;
        enth_lst.omega = omega;
        enth_lst.CTS = CTS;

        % program termination
        if any(isnan(u(:))) || any(isnan(E(:)))
            break
        end

    end

    % print results
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,is_active)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,is_active)))

    [SMB, smb_driver] = run_smb( ...
        hS, is_active, smb_driver, arrayTime(iTimeStep), p);

    % Store the state at the current output time before projection.
    At_u(:,:,iTimeStep) = u;
    At_hS(iTimeStep, :) = hS;
    At_hB(iTimeStep, :) = hB;
    At_H(iTimeStep, :) = H;
    At_SMB(iTimeStep, :) = SMB;

    %% ICE THICKNESS EVOLUTION (adaptive CFL subcycling)
    if p.is_thk_evolv && iTimeStep < numTimeStep
        [dt_CFL, umax] = calc_thk_CFL_step(u_s, p, dx);
        t_local = 0;
        nsub = 0;
        nsub_max = ceil(dt_u/p.dt_H_min) + 1;

        while t_local < dt_u - 1e-12
            remaining = dt_u - t_local;
            dt_H = min([dt_CFL, p.dt_H_max, remaining]);

            nsub = nsub + 1;
            if nsub > nsub_max
                error('main_polyT:ExcessiveThicknessSubsteps', ...
                    ['Thickness evolution exceeded %d substeps over %.6g yr. ' ...
                     'The geometry or velocity solution is likely unstable.'], ...
                    nsub_max, dt_u)
            end

            t_next = t_local + dt_H;
            if ~isfinite(t_next) || t_next <= t_local
                error('main_polyT:ThicknessTimeStalled', ...
                    ['Thickness time failed to advance: t = %.16g yr and ' ...
                     'dt = %.6g yr.'], t_local, dt_H)
            end

            Hnp1 = get_evolution_continuity(H, u, SMB, dt_H);

            i_term_old = find(is_active, 1, 'last');
            while ~isempty(i_term_old) && Hnp1(i_term_old) <= p.Hmin
                is_active(i_term_old) = false;
                i_term_old = find(is_active, 1, 'last');
            end

            Hnp1(Hnp1 < p.Hmin) = p.Hmin;
            H = Hnp1;
            hS = hB + H;
            t_local = t_next;
            SMB(~is_active) = 0;

            fprintf('  thickness substep dt = %.6g yr, umax = %.6g m/yr\n', ...
                dt_H, umax)

            if ~any(is_active)
                fprintf('Glacier became extinct during thickness evolution.\n')
                break
            end

            % Thermodynamics remains operator-split over dt_u, but the
            % geometry-dependent velocity and CFL limit are refreshed.
            if t_local < dt_u - 1e-12
                set_staggered_grid();
                [u_s, u, visc_s, visc, strainHeat] = ...
                    solve_u_smedt(visc_s, visc, AGlen_s, p);
                [dt_CFL, umax] = calc_thk_CFL_step(u_s, p, dx);
            end
        end

        if ~any(is_active)
            break
        end
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

