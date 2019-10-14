clc
clearvars
clearvars -global

tic
global SPY N M Ms iter_max iter_u u_s_lst
set_ice_parameters;
set_ice_geometry;
set_staggered_grid();

% Initialization
% AGlen_s = 1e-16*ones(N, Ms); % [Pa-3 a-1]
visc_s = 1e13/SPY*ones(N, Ms); % [Pa a]
visc = 1e13/SPY*ones(N, M); % [Pa a]

iter_u = 0;
while 1
    iter_u = iter_u + 1;
    fprintf('iter_u: %d\n', iter_u)
    
    u_s = solver_u(visc_s, visc);
    
    % convergence
    if iter_u == 1
        gamma = 1;
        gamma_min = 0.8;
        rho_t = 0.9;
        tol = 1e-7;
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
    
    u = staggerX2main(u_s);   
    [visc_s, visc, ~] = get_ice_viscosity(u_s, u);
end

toc
fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))

plot_u_field_us