clc
clearvars
clearvars -global

addpath('..')

tic
global N M Ms iter_u u_s_lst beta2_s usurf_s_obs

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
p.type_BBC = 'LFlaw'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT';

%% GLACIER GEOMETRY
%
%
set_ice_geometry('../geo_inputs/geo_arolla', p);
obs = load('obs_forward_arolla_1e4_flowband0.mat');
beta2_obs = obs.beta2_obs;
usurf_obs = obs.usurf_obs;
ub_obs = obs.ub_obs;

usurf_s_obs = main2staggerX(usurf_obs);
beta2_s_obs = main2staggerX(beta2_obs);
ub_s_obs = main2staggerX(ub_obs);

%% ROBIN INVERSION
%
%
% WHILE LOOP with updated para_mu
para_mu = 3;
para_mu0 = para_mu;
is_covg_Robin = 0;
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

% initial guess beta2
beta2 = 1e4*ones(1,M);
beta2_s = main2staggerX(beta2);

% WHILE LOOP with specified para_mu
cvg_cost_arr = [];
iter_robin = 0;
while 1
    iter_robin = iter_robin + 1;
%     fprintf('para_mu: %3.2f\n', para_mu)
    
    % switch the surface boundary condition
    if mod(iter_robin,2)~=0
        p.type_SBC = 'neum';
    else
        p.type_SBC = 'diri';
    end
    fprintf('iter_robin: %d, para_mu: %3.2f, type_SBC: %s\n', iter_robin, para_mu, p.type_SBC)
    
    set_staggered_grid();
    
    %% velocity solver
    %
    %
    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        
        [u_s] = solver_u_core(visc_s, visc, AGlen_s, p);
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
    end
    % calculate the vertical velocity
    [w_vs, w] = get_ice_w(u_s, u);
    
    fprintf('usurf_mean: %3.1f, usurf_min: %3.8f, usurf_max: %3.1f \n',...
        [mean(u(end,:)), min(u(end,:)), max(u(end,:)), ])
    %% Robin inversion
    %
    %
    if strcmp(p.type_SBC, 'neum')
        ub_neum = u(1,:);
        uN = u;
    else
        ub_diri = u(1,:);
        uD = u;
        
        beta2 = beta2.*((abs(ub_neum)+0.1)./(abs(ub_diri)+0.1)).^para_mu;
        
        beta2(beta2>1e6) = 1e6; % upper bound 1e6 Pa a m-1 (Schafer et al., 2012)
        beta2(beta2<0) = 1; % lower bounda 1 Pa a m-1 (Schafer et al., 2012)
        beta2_s = main2staggerX(beta2);
        
        cost_robin = sum(abs(uN(end,:) - uD(end,:)));
        if iter_robin>2            
            cvg_cost = abs(cost_robin-cost_robin_lst) / cost_robin;
            cvg_cost_arr = [cvg_cost_arr; cvg_cost];
            if cvg_cost<5e-2 % /cost_robin_lst
                is_covg_Robin = 1;
                disp('Robin inversion converged!')
                break
            end
        end
        cost_robin_lst = cost_robin;
    end
    
    if iter_robin>30
        break
    end
    
    para_mu = para_mu - 0.02;
end

toc

% save result_inv xi xi_s u u_s beta2 beta2_s cost_arr para_mu
% save result_inv_ERG

%% Figure
%
global xi
xi_s = main2staggerX(xi);
figure
subplot(3,1,1)
hold on
plot(xi_s/1e3, usurf_s_obs, 'k-', 'linewidth', 3);
plot(xi_s/1e3, u_s(end,:), 'r-', 'linewidth', 0.5,...
    'Marker', '.', 'MarkerSize', 10);
hold off
title(strcat('\nu_{init}=', num2str(para_mu0), ' ,\nu_{end}=',num2str(para_mu)))
xlabel('Horizontal distance (km)')
ylabel('u_s (m a^{-1})')
hl = legend('u_s observed', 'u_s modeled');
set(hl, 'location', 'NorthEastOutside')
box on

subplot(3,1,2)
hold on
plot(xi_s/1e3, beta2_s_obs, 'k-', 'linewidth', 3);
plot(xi_s/1e3, beta2_s, 'r-', 'linewidth', 0.5,...
    'Marker', '.', 'MarkerSize', 10);
ylim([0,max(beta2_s_obs)+1e3])
hold off

xlabel('Horizontal distance (km)')
ylabel('\beta^2 (Pa m^{-1} a)')
box on

hl = legend('\beta^2 observed', '\beta^2 modeled');
set(hl, 'location', 'NorthEastOutside')

subplot(3,1,3)
hold on
plot(xi_s/1e3, ub_s_obs, 'k-', 'linewidth', 3);
plot(xi_s/1e3, u_s(1,:), 'r-', 'linewidth', 0.5,...
    'Marker', '.', 'MarkerSize', 10);
hold off
xlabel('Horizontal distance (km)')
ylabel('u_b (m a^{-1})')
hl = legend('u_b observed', 'u_b modeled');
set(hl, 'location', 'NorthEastOutside')
box on

set(gcf, 'position', [300, 300, 560, 420])