clc
clearvars
clearvars -global

tic
global SPY N Ms iter_max iter_u u_s_lst
set_ice_parameters;
set_ice_geometry;
set_staggered_grid();

% Initialization
AGlen_s = zeros(N, Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N, Ms) + 1e13/SPY; % [Pa a]

iter_u = 0;
while 1
    iter_u = iter_u + 1;
    fprintf('iter_u: %d\n', iter_u)
    
    u_s = solver_u_classic(visc_s);
    
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
    [visc_s, visc, ~] = get_ice_viscosity_classic(u_s, AGlen_s);
end

[tauxz, tauxz_s] = calc_tauxz(u_s, visc_s);
[Pdiff, Pdiff_s] = calc_Pdiff(u_s, visc_s);

toc
fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))

plot_u_field_us
figure
plot(xi/1000, tauxz)