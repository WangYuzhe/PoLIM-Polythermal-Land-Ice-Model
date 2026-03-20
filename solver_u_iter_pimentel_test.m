function [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_pimentel_test(visc_s,visc,AGlen_s,para)
% Wrapped Picard iteration for velocity solver
% Based on the adaptive underrelaxation subspace solver proposed by Smedt et al. (2010)
% Created on 2023/12/1
% Updated on 2025/11/1

global iter_u u_s_lst

iter_u = 0;

while iter_u<=para.iter_max
    iter_u = iter_u + 1;
    fprintf('iter_u: %d\n', iter_u)

    % --- Core velocity solver ---
    u_s = solver_u_core_test(visc_s, visc, AGlen_s, para);

    % --- Begin Picard iteration ---
    if iter_u == 1
        gamma_cvg = 1.0;
        gamma_cvg_min = 0.2; % results are not sensitive to this variable
        rho_cvg = 0.9;
        tol_cvg = 1e-5;
    elseif iter_u > 1 && iter_u < para.iter_max
        u_s_now = u_s;
        delta_cvg = norm(u_s_now - u_s_lst);
        gamma_cvg_min = rho_cvg*gamma_cvg_min;
        alpha_cvg = -log((gamma_cvg*rho_cvg - gamma_cvg_min)/...
            (1 - gamma_cvg_min))/(delta_cvg - tol_cvg);
        if delta_cvg > tol_cvg
            gamma_cvg = gamma_cvg_min + (1-gamma_cvg_min)*...
                exp(-alpha_cvg*(delta_cvg-tol_cvg));
        else
            gamma_cvg = 1.0;
        end
        u_s = u_s_lst + gamma_cvg*(u_s_now - u_s_lst);

        if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now) < 5e-4
            break
        end
    else
        break
    end

    u_s_lst = u_s;

    % --- Update viscosity and strain heating ---
    u = staggerX2main(u_s);
    [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, para);
end
end