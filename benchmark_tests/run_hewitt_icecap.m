clc
clearvars
clearvars -global

addpath('..')

tic
global N M hB H iTimeStep zeta dhSdx enth_lst

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params_hewitt_icecap();
SPY = p.SPY;
p.Hmin = 0;

p.type_enth = 'MEGM';
p.is_auto_enth_BBC = 0; % At bed: T = Tm
p.has_Greve_drainage = 1;
p.layers = 201;

%% GLACIER GEOMETRY
%
%
set_ice_geometry('../geo_inputs/geo_hewitt_icecap', p);

%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 1000; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, endTime);

%% INITIALIZATION
%
%
%-----------------solve the velocity field-----------------
set_staggered_grid();

u = zeros(N,M);
strainHeat = zeros(N,M);
for i = 1:M
    for j = 1:N
        u(j,i) = -2*p.AGlen*(p.rhoi*p.g*dhSdx(i)).^3/4*(H(i).^4 - ((1-zeta(j))*H(i)-hB(i)).^4);
        strainHeat(j,i) = 2*p.AGlen*(p.rhoi*p.g*dhSdx(i)).^4*((1-zeta(j)).*H(i)-hB(i)).^4;
    end
end
u_s = main2staggerX(u);
[w_vs, w] = get_ice_w(u_s, u);

% convert [m s-1] to [m a-1]
u = u*SPY;
u_s = u_s*SPY;
w = w*SPY;
w_vs = w_vs*SPY;
strainHeat = strainHeat*SPY;
frictionHeat = zeros(1,M);
%-----------------solve the velocity field-----------------

%  Related to enthalpy solver
enth_covg = 0;
enth_nan = 0;
pdtj_storage = [];

% SBC and IC for the enthalpy solver
Ts = -1;
Esbc = p.Cp*(Ts + 273.15 - p.Tref)*ones(1,M);
Eini = p.Cp*(262.15 - p.Tref)*ones(N,M);

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    % calculate the enthalpy
    if strcmpi(p.type_enth,'SEGM')
        % [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP_drain_Greve, qw_TEMP_diffu] = ...
        %     solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
        
        [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP_drain_greve, qw_TEMP_diffu] = ...
            solver_enthalpy_SEGM_blatter(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
    elseif strcmpi(p.type_enth,'MEGM')
        [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP, qw_TEMP_darcy] =...
            solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
%         
%         [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP, qw_TEMP_darcy] =...
%             solver_enthalpy_MEGM_blatter(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
    elseif strcmpi(p.type_enth,'isoT')
        % pass
    end
    
    % converge
    %     if iTimeStep>1
    %         Elst = At_E(:,:,iTimeStep-1);
    %         pdtj = sumsqr(E-Elst)/sumsqr(Elst);
    %         pdtj_storage = [pdtj_storage; pdtj];
    %         if pdtj<2.1e-8
    %             enth_covg = 1;
    %             break
    %         end
    %     end
    
    %     if iTimeStep>1
    %         epsilon_enth = sumsqr(E-E_lst)/sumsqr(E_lst);
    %         if epsilon_enth<1e-8
    %             break
    %         end
    %     end
    %     E_lst = E;
    
    %%
    %
    % STORAGE RESULTS
    enth_lst.E = E;
    enth_lst.omega = omega;
    enth_lst.Kappa_vs = Kappa_vs;
    enth_lst.thk_w = thk_w;
    enth_lst.thk_TEMP = thk_TEMP;
    enth_lst.is_TEMP = is_TEMP;
    
    % logic value for basal thermal state
    
    if strcmpi(p.type_enth,'MEGM')
        enth_lst.qw_TEMP_darcy = qw_TEMP_darcy;
    end
end
toc

plot_enthalpy_Hewitt_icecap
