clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms hS hB H xi iTimeStep

%% PHYSICAL CONSTANTS AND PARAMETERS
%
p = params();
SPY = p.SPY;

p.Hmin = 0.1;
p.layers = 15;
p.type_valley = 'rect';

p.is_flowband = 0;
p.is_thk_evolv = 1;
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
endTime = 50; % end time of model run [a]
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
At_smb = zeros(numTimeStep,M);

%% MAIN
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d, year: %d \n', iTimeStep, arrayTime(iTimeStep))

    set_staggered_grid();

    dx = mean(diff(xi)); % assuming x is your horizontal grid

    [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_smedt(visc_s, visc, AGlen_s, p);
    %% ICE THICKNESS EVOLUTION (adaptive CFL time stepping)
    if p.is_thk_evolv
        t_local = 0.0; % local time within one velocity step

        while t_local < dt_u - 1e-12
            % --- Compute SMB (can be updated every substep) ---
            grad_smb = 1.5e-2;
            smb = calc_smb_gradient(hS, median(hS), grad_smb);
            smb(H <= p.Hmin) = 0;

            % --- CFL-based adaptive time step ---
            umax = max(abs(u(:)));

            if umax > 0
                dt_CFL = p.CFL * dx / umax;
            else
                dt_CFL = p.dt_H_max; % no flow → large step allowed
            end

            % Apply bounds
            dt_H = min([dt_CFL, p.dt_H_max, dt_u - t_local]);
            dt_H = max(dt_H, p.dt_H_min);

            % --- Thickness update ---
            Hnp1 = get_evolution_continuity(H, u, smb, dt_H);

            % Enforce minimum thickness
            Hnp1(Hnp1<p.Hmin) = p.Hmin;

            % Update state
            H = Hnp1;
            hS = hB + H;

            % Advance local time
            t_local = t_local + dt_H;

            fprintf('  substep dt = %.4f yr, umax = %.2f m/yr\n', dt_H, umax);
        end

        % Store annual SMB (optional: last substep or mean)
        At_smb(iTimeStep, :) = smb;
    end

    if p.is_surf_relax
        Hn = H;
        hSn = hS;
        SMB = zeros(1,M);
        SMB(Hn<p.Hmin) = 0;
        At_smb(iTimeStep, :) = SMB;

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

