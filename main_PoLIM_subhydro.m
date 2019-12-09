clc
clearvars
clearvars -global

global SPD SPY N M Ms iter_max iter_u u_s_lst iTimeStep
set_ice_parameters;
set_ice_geometry();

%% Time setting
dt = 1; % [a]
endTime = 1; % [a]
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

dt_hydro = 1*SPD;
endTime_hydro = 1*SPY;
[arrayTime_hydro, numTimeStep_hydro] = set_time_step(dt_hydro, endTime_hydro);

%% Initialization
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

%  Related to velocity solver
At_u = zeros(N,M,numTimeStep);
At_w = zeros(N,M,numTimeStep);

%  Related to glacier evolution
At_hS = zeros(numTimeStep,M);
At_hB = zeros(numTimeStep,M);

%  Related to subglacier hydrology
At_Hw = zeros(numTimeStep, M);
At_dHw3dx = zeros(numTimeStep,M);
At_phi = zeros(numTimeStep,M);
At_N = zeros(numTimeStep, M);
At_qw = zeros(numTimeStep, M);
iter_subhydro = 0;
iter_subhydro_threshold = 1;
%% MAIN
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    set_staggered_grid();
    
    %-----------subglacier hydrology-----------
    if iTimeStep == 1
        load result_A1.mat
        phi_nm1 = phi_A1;
        Hw_nm1 = Hw_A1; % [m]
        dHw3dx_nm1 = dHw3dx_A1;
    else
        Hw_nm1 = At_Hw(iTimeStep-1,:);
        dHw3dx_nm1 = At_dHw3dx(iTimeStep-1,:);
        phi_nm1 = At_phi(iTimeStep-1,:);
    end
    
    [Neff, phi, qw, Hw, dHw3dx] = subHydro_storage_Hoffman(dt,...
        basalMelt, mouInput, seasonRunoff(iTimeStep,:), Hw_nm1, dHw3dx_nm1, phi_nm1, ev, ub, visc_ice_b);
    
    At_Hw(iTimeStep, :) = Hw;
    At_dHw3dx(iTimeStep, :) = dHw3dx;
    At_phi(iTimeStep, :) = phi;
    At_N(iTimeStep, :) = Neff;
    At_qw(iTimeStep, :) = qw;
    
    Neff_s = main2staggerX(Neff);
    %-----------subglacier hydrology-----------
    
    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        fprintf('iter_u: %d\n', iter_u)
        
        [u_s] = solver_u_subHydro(visc_s, visc, AGlen_s, Neff_s);
        %-------------------------Picard iteration-------------------------
        rho_t = 0.9;
        tol = 1e-7;
        if iter_u == 1
            gamma = 1;
            gamma_min = 0.8;
        elseif iter_u > 1 && iter_u < iter_max
            u_s_now = u_s;
            sigma = norm(u_s_now - u_s_lst);
            gamma_min = rho_t*gamma_min;
            alpha = -log((gamma*rho_t-gamma_min)/(1-gamma_min))/(sigma-tol);
            if sigma > tol
                gamma = gamma_min + (1-gamma_min)*exp(-alpha*(sigma-tol));
            else
                gamma = 1;
            end
            u_s = u_s_lst + gamma*(u_s_now - u_s_lst);
            
            if sigma/norm(u_s_now) < 1e-2
                break
            end
        else
            break
        end
        u_s_lst = u_s;
        %-------------------------Picard iteration-------------------------
        u = staggerX2main(u_s);
        
        % calculate the vertical velocity field
        [w_vs, w] = get_ice_w(u_s, u);
        
        % calculate the strain heat field
        [visc_s, visc, ~] = get_ice_viscosity(u_s, u, AGlen_s);
    end
    
    At_u(:,:,iTimeStep) = u;
    At_w(:,:,iTimeStep) = w;
    
    
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
end

% plot_enthalpy_field
plot_uField_uSurf