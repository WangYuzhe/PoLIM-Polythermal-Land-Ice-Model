clc
clearvars
clearvars -global

addpath('.')
addpath('.\mini_tools')

tic
global N M Ms dx hS hB H iTimeStep is_active

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
is_active = H > p.Hmin;

%% TIME SETTING
%
dt_u = 1.0; % time step for velocity solver [a]
startTime = 1; % start time of model run [a]
endTime = 1; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, startTime, endTime);

p.CFL = 0.4;
p.dt_H_max = dt_u;
p.dt_H_min = 0.001;

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

    if ~any(is_active)
        fprintf('Glacier is extinct; stopping the simulation.\n')
        break
    end

    set_staggered_grid();

    % [u_s, u, visc_s, visc, strainHeat] = solve_u_smedt(visc_s, visc, AGlen_s, p);
    [u_s, u, visc_s, visc, strainHeat] = solve_u_pimentel(visc_s, visc, AGlen_s, p);
    
    zela = median(hS(is_active));
    grad_smb = 5.0e-2;
    m_b = calc_smb_gradient(hS, zela, grad_smb);
    m_b(~is_active) = 0;

    % Store the state at arrayTime(iTimeStep), before projecting it to the
    % next output time.
    At_u(:,:,iTimeStep) = u;
    At_hS(iTimeStep, :) = hS;
    At_hB(iTimeStep, :) = hB;
    At_H(iTimeStep, :) = H;
    At_SMB(iTimeStep, :) = m_b;

    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,is_active)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,is_active)))

    if any(isnan(u(:))) || any(H > 1000)
        warning('main_isoT_forward:UnstableState', ...
            'The model became unstable at time step %d.', iTimeStep)
        break
    end

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
                error('main_isoT_forward:ExcessiveThicknessSubsteps', ...
                    ['Thickness evolution exceeded %d substeps over %.6g yr. ' ...
                     'The geometry or velocity solution is likely unstable.'], ...
                    nsub_max, dt_u)
            end

            t_next = t_local + dt_H;
            if ~isfinite(t_next) || t_next <= t_local
                error('main_isoT_forward:ThicknessTimeStalled', ...
                    ['Thickness time failed to advance: t = %.16g yr and ' ...
                     'dt = %.6g yr.'], t_local, dt_H)
            end

            Hnp1 = get_evolution_continuity(H, u, m_b, dt_H);

            % Retreat only from the downstream end; an interior thin cell
            % must not split the flowline.
            i_term_old = find(is_active, 1, 'last');
            while ~isempty(i_term_old) && Hnp1(i_term_old) <= p.Hmin
                is_active(i_term_old) = false;
                i_term_old = find(is_active, 1, 'last');
            end

            Hnp1(Hnp1 < p.Hmin) = p.Hmin;
            H = Hnp1;
            hS = hB + H;
            t_local = t_next;
            m_b(~is_active) = 0;

            fprintf('  thickness substep dt = %.6g yr, umax = %.6g m/yr\n', ...
                dt_H, umax)

            if ~any(is_active)
                fprintf('Glacier became extinct during thickness evolution.\n')
                break
            end

            % Keep geometry, velocity, and the CFL limit mutually
            % consistent throughout the projection interval.
            if t_local < dt_u - 1e-12
                set_staggered_grid();
                [u_s, u, visc_s, visc, strainHeat] = ...
                    solve_u_pimentel(visc_s, visc, AGlen_s, p);
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

% plot the results
plot_uField_uSurf

