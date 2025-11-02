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

p.is_flowband = 1;
p.is_thk_evolv = 0;
p.is_surf_relax = 1;

p.type_BBC = 'CFlaw_isoT'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'
p.lambda_max = 4; % wavelength of the dominant bedrock bumps [m]
p.m_max = 0.3; % maximum slope of the dominant bedrock bumps []; ref: 0.5

%% GLACIER GEOMETRY
%

file_geo = './geo_inputs/geo_lhg12_glate.mat';
set_ice_geometry(file_geo, p);

glac = load(file_geo);
geo = glac.geo;

%% TIME SETTING
%
dt_u = 0.1; % time step for velocity solver [a]
start_yr = 1;
end_yr = 3;
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

    % calculate the horizontal velocity    
    [u_s, u, visc_s, visc, strainHeat] = solver_u_iter(visc_s, visc, AGlen_s, p);
    
    % calculate the vertical velocity
    [w_vs, w] = get_ice_w(u_s, u);
    %% ICE THICKNESS EVOLUTION
    %
    %
    if p.is_thk_evolv
        Hn = H;
        hSn = hS;

        zela = 300; % median(hS(1:40))
        grad_smb = 4.0e-2;
        smb = calc_smb_gradient(hSn, zela, grad_smb);

        smb(Hn<p.Hmin) = 0;
        At_SMB(iTimeStep, :) = smb;

        Hnp1 = get_evolution_continuity(Hn, u, smb, dt_u);
        % Hnp1 = get_evolution_continuity_chatgpt(Hn, u, smb, dt_u);

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

if p.is_surf_relax
    geo.hS = hS;
    geo.H = H;
    save('./geo_inputs/geo_lhg12_glate_relax', 'geo')
end

% plot the results
figure
plot_uField_uSurf

