clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms hS hB H dzeta iter_u u_s_lst iTimeStep...
    At_E At_T At_omega At_CTS At_Kappa_vs...
    At_isTEMP At_thk_TEMP At_thk_w At_qw_TEMP_darcy

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
p.is_thk_evolv = 0;
p.is_surf_relax = 0;
p.is_subHyro = 0;
p.type_BBC = 'LFlaw_E2'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'; 'CFlaw_polyT'

p.type_Arrhenius = 'cuffey';
p.type_enth_solver = 'isoT';
p.is_auto_enth_BBC = 1;
p.type_enth_BBC = 2;
p.is_duval = 1;
p.moi_UB = 0.01;
p.k0 = 1e-12;

p.Qgeo = 0.06; % 0.06, ref
p.kflot = 0.1;
%% GLACIER GEOMETRY
%
set_ice_geometry('../geo_inputs/geo_arolla_resample', p);

%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 1; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, endTime);

%% INITIALIZATION
%
%
% Related to Ney-Glen law
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]

visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

% Related to velocity solver
At_u = zeros(N,M,numTimeStep);
At_w = zeros(N,M,numTimeStep);

% Related to enthalpy solver
At_E = zeros(N,M,numTimeStep);
At_T = zeros(N,M,numTimeStep);
At_omega = zeros(N,M,numTimeStep);
At_CTS = zeros(numTimeStep,M);
At_Kappa_vs = zeros(N-1,M,numTimeStep);
At_isTEMP = zeros(numTimeStep,M);
At_thk_w = zeros(numTimeStep,M); % water thickness
At_thk_TEMP = zeros(numTimeStep,M);
At_qw_TEMP_darcy = zeros(N-1,M,numTimeStep);
At_m_basal = zeros(numTimeStep,M); % basal melt rate

% Related to glacier evolution
At_hS = zeros(numTimeStep,M);
At_hB = zeros(numTimeStep,M);
At_H = zeros(numTimeStep,M);
At_SMB = zeros(numTimeStep,M);

% SBC and IC for the enthalpy solver
if ~strcmpi(p.type_enth_solver,'isoT')
    % Thermal surface boundary condition
    Esbc = set_thermalSBC(p);
    
    % Initial enthalpy field
    Eini = get_initial_enthalpy(Esbc, p);
end

%% MAIN
%
%
pdtj_storage = [];
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    set_staggered_grid();
    
    if ~strcmpi(p.type_enth_solver,'isoT')
        % Thermal surface boundary condition
        Esbc = set_thermalSBC(p);
    end
    
    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        fprintf('iter_u: %d\n', iter_u)
        
        [u_s] = solver_u(visc_s, visc, AGlen_s, p);
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
        
        if iter_u>=iter_max
            break
        end
        
        u_s_lst = u_s;
        %-----------------------End Picard iteration-----------------------
        u = staggerX2main(u_s);
        
        % calculate the vertical velocity
        [w_vs, w] = get_ice_w(u_s, u);
        
        % calculate the strain heat
        [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, p);
        [tauxz, ~] = calc_tauxz(u_s, visc_s);        
        frictionHeat = tauxz.*u(1,:)/p.SPY; % [W m-2]
        
        % calculate the enthalpy
        if strcmpi(p.type_enth_solver,'SEGM')
            [E, T, omega, Kappa_vs, CTS, thk_TEMP, m_basal, qw_greveDrain, qw_TEMP_diffu] = ...
                solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
        elseif strcmpi(p.type_enth_solver,'MEGM')
            [E, T, omega, Kappa_vs, CTS, thk_TEMP, m_basal, qw_TEMP, qw_TEMP_darcy] =...
                solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
        elseif strcmpi(p.type_enth_solver,'isoT')
            % pass
        end
        
        % calculate the AGlen
        if ~strcmpi(p.type_enth_solver,'isoT')
            AGlen_s = get_AGlen(T, omega, CTS, p);
        end
    end    
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
        
        Hnp1 = get_evolution_continuity(Hn, u, SMB, dt_u);
        Hnp1(Hnp1<p.Hmin) = p.Hmin;
        hS = hB + Hnp1;
        H = Hnp1;
    end
    
    if p.is_surf_relax
        Hn = H;
        %         SMB = 5e-3*(hS - 3150); % 3940 is for Arolla; 4000~4100 for incline
        hSn = hS;
        zref = 2500;
        %         SMB = calc_SMB(zref, Tma, hSn);
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
    % STORAGE RESULTS
    At_u(:,:,iTimeStep) = u;
    At_w(:,:,iTimeStep) = w;
    At_hS(iTimeStep, :) = hS;
    At_hB(iTimeStep, :) = hB;
    At_H(iTimeStep, :) = H;
    
    if ~strcmpi(p.type_enth_solver,'isoT')
        At_E(:,:,iTimeStep) = E;
        At_T(:,:,iTimeStep) = T;
        At_omega(:,:,iTimeStep) = omega;
        At_Kappa_vs(:,:,iTimeStep) = Kappa_vs;
        At_CTS(iTimeStep,:) = CTS;
        At_thk_TEMP(iTimeStep,:) = thk_TEMP;
        At_m_basal(iTimeStep,:) = m_basal;
        
        % calculate the basal water layer thickness
        if iTimeStep == 1
            At_thk_w(iTimeStep,:) = zeros(1,M) + dt_u*m_basal;
        else
            At_thk_w(iTimeStep,:) = At_thk_w(iTimeStep-1,:) + dt_u*m_basal;
        end
        
        % logic value for basal thermal state
        logic1 = (At_thk_w(iTimeStep,:)>0 & At_isTEMP(iTimeStep,:)==0);
        At_isTEMP(iTimeStep,:) = At_isTEMP(iTimeStep,:) | logic1;
    end
    
    if strcmpi(p.type_enth_solver,'MEGM')
        At_qw_TEMP_darcy(:,:,iTimeStep) = qw_TEMP_darcy;
    end
    
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
    
end
toc

global xi
normalized_xi = (xi - xi(1))./(xi(end) - xi(1));

figure
set(gcf, 'position', [146,332,1000,350]);
subplot(1,2,1)
load('ismip_hom_us_mean_std.mat')
plot_benchmark.plot_shadedErrorBar(NFS_e001, FS_e001, 'Velocity (m a^{-1})')
hold on
plot(normalized_xi, u(end,:), 'k--', 'linewidth', 1)

subplot(1,2,2)
load('ismip_hom_expE_taub_mean_std.mat')
plot_benchmark.plot_shadedErrorBar(NFS_e001, FS_e001, 'Basal shear stress (kPa)')
hold on
plot(normalized_xi, tauxz/1e3, 'k--', 'linewidth', 1)