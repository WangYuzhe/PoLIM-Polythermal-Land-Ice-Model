clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms iter_u u_s_lst iTimeStep

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params();
SPY = p.SPY;
iter_max = p.iter_max;
p.Hmin = 0.1;
p.layers = 51;
p.type_valley = 'rect';
p.type_BBC = 'zero'; %'LFlaw';
p.is_flowband = 0;

%%
%
%
set_staggered_grid_period();
%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 1; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, 1, endTime);

%% INITIALIZATION
%
%
% Related to Ney-Glen law
AGlen_s = 1e-16*ones(N,Ms); % [Pa-3 a-1]
visc_s = 1e13/SPY*ones(N,Ms); % [Pa a]
visc = 1e13/SPY*zeros(N,M); % [Pa a]

% Related to velocity solver
At_u = zeros(N,M,numTimeStep);
At_w = zeros(N,M,numTimeStep);

%% MAIN
%
%

for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)

    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        fprintf('iter_u: %d\n', iter_u)
        
        [u_s] = solver_u_period(visc_s, visc, AGlen_s, p);
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
        u = (u_s(:,1:end-1) + u_s(:,2:end))/2;       
        
        % calculate the vertical velocity
        [w_vs, w] = get_ice_w(u_s, u);
        
        % calculate the strain heat
        [visc_s, visc, strain_heat] = get_ice_viscosity_period(u_s, u, AGlen_s, p);
    end
end

toc

% plot_uField_uSurf

global xi
normalized_xi = (xi - xi(1))./(xi(end) - xi(1));

figure
% set(gcf, 'position', [146,332,1000,350]);
load('ismip_hom_us_mean_std.mat')
plot_benchmark.plot_shadedErrorBar(NFS_b160, FS_b160, 'Velocity (m a^{-1})')
hold on
plot(normalized_xi, u(end,:), 'k--', 'linewidth', 1)