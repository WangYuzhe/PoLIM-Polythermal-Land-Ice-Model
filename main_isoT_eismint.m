clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms hS hB H iter_u u_s_lst iTimeStep

%% PHYSICAL CONSTANTS AND PARAMETERS
%
p = params();
SPY = p.SPY;

p.Hmin = 0.1;
p.layers = 25;
p.type_valley = 'rect';

p.is_flowband = 0;
p.is_thk_evolv = 1;

%% GLACIER GEOMETRY
%
set_ice_geometry('./EISMINT/geo_eismint.mat', p);
% set_ice_geometry('./geo_inputs/geo_arolla.mat', p);

%% TIME SETTING
%
dt_u = 100; % time step for velocity solver [a]
start_yr = 1;
end_yr = 10e3;
[arrayTime, numTimeStep] = set_time_step(dt_u,start_yr,end_yr);

%% INITIALIZATION
%
% Related to Ney-Glen law
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]

visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

% Related to velocity solver
At_u = zeros(N,M,numTimeStep);
At_w = zeros(N,M,numTimeStep);

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

    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        fprintf('iter_u: %d\n', iter_u)

        [u_s] = concise_solver_u(visc_s, visc, AGlen_s, p);

        % [u_s] = solver_u(visc_s, visc, AGlen_s, p);
        %----------------------Begin Picard iteration----------------------
        if iter_u>2
            u_s_now = u_s;
            Cs = u_s_now - u_s_lst;
            Sita = acos(Cs'*C/(sumsqr(Cs)*sumsqr(C)));
            if isequal(Sita<=pi/8, ones(Ms,Ms))
                mu1 = 2.5;
            elseif isequal(Sita>pi/8,ones(Ms,Ms)) && isequal(Sita<19*pi/20, ones(Ms,Ms))
                mu1 = 1;
            elseif isequal(Sita>=19*pi/20, ones(Ms,Ms))
                mu1 = 0.5;
            end
            u_s = u_s_lst + mu1*Cs;

            if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now)<1e-4
                break
            end
        end

        if iter_u>1
            C = u_s - u_s_lst;
        end

        if iter_u>=p.iter_max
            break
        end

        u_s_lst = u_s;
        %-----------------------End Picard iteration-----------------------
        u = staggerX2main(u_s);

        % calculate the vertical velocity
        [w_vs, w] = get_ice_w(u_s, u);

        % calculate the strain heat
        [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, p);
    end
    %% ICE THICKNESS EVOLUTION
    %
    %
    if p.is_thk_evolv
        Hn = H;
        hSn = hS;

        smb = 0.3*ones(1,M);

        smb(Hn<p.Hmin) = 0;
        At_SMB(iTimeStep, :) = smb;

        % Hnp1 = get_evolution_continuity(Hn, u, smb, dt_u);
        Hnp1 = get_evolution_continuity_chatgpt(Hn, u, smb, dt_u);
        % Hnp1 = get_evolution_ai(Hn, u, smb, dt_u);

        Hnp1(Hnp1<=p.Hmin) = p.Hmin;
        hS = hB + Hnp1;
        H = Hnp1;
    end
    %%
    %
    % RESULTS
    At_u(:,:,iTimeStep) = u;
    At_w(:,:,iTimeStep) = w;
    At_hS(iTimeStep, :) = hS;
    At_hB(iTimeStep, :) = hB;
    At_H(iTimeStep, :) = H;

    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))

    %%
    if sum(isnan(u(:)))
        break
    end

    if sum(H<=p.Hmin)==M
        break
    end

end
toc

% plot the results
plot_uField_uSurf

figure;
plot(xi, At_hS(1, :), 'k-')
hold on
plot(xi, At_hS(end, :), 'r-')

