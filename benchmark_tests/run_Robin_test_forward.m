clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms iter_u u_s_lst iTimeStep xi beta2_s

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params();
SPY = p.SPY;
iter_max = p.iter_max;
p.Hmin = 0.1;
p.layers = 31;
p.type_valley = 'rect';

p.is_flowband = 0;
p.type_BBC = 'LFlaw'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT'; 'CFlaw_polyT'
p.type_SBC = 'neum';

%% GLACIER GEOMETRY
%
set_ice_geometry('../geo_inputs/geo_arolla', p);

%% Prescribed basal friction parameter
%
beta2 = 1e4 + 1e4*sin(2*pi/xi(end)*xi);
beta2_s = main2staggerX(beta2);
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
        
        % calculate the strain heat
        [visc_s, visc, ~] = get_ice_viscosity(u_s, u, AGlen_s, p);
        
        disp(norm(u, 'fro'))
    end
    
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
end

% calculate the vertical velocity
[w_vs, w] = get_ice_w(u_s, u);
    
    

%% Save results
%
usurf_obs = u(end,:);
beta2_obs = beta2;
ub_obs = u(1,:);
save('obs_forward_arolla_1e4_flowband0.mat', 'usurf_obs', 'ub_obs', 'beta2_obs')

%% Figure
%
plot_uField_uSurf

figure
set(gcf, 'position', [178,317,1269,242]);
subplot(1,2,1)
plot(xi/1e3, u(end,:), 'k-', 'linewidth', 1)
xlabel('Horizontal distance (km)')
ylabel('Surface velocity (m a^{-1})')

subplot(1,2,2)
plot(xi/1e3, beta2, 'k-', 'linewidth', 1)
xlabel('Horizontal distance (km)')
ylabel('\beta (Pa a m^{-1})')