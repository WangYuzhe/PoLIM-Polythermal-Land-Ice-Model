function [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_chatgpt(visc_s,visc,AGlen_s,para)
% Wrapped Picard iteration for velocity solver
% Based on the relaxed Picard algorithm proposed by De Smedt et al. (2010)
% Created on 2023/12/1
% Updated on 2025/11/1

global iter_u u_s_lst

iter_u = 0;

while iter_u<=para.iter_max
    iter_u = iter_u + 1;
    fprintf('iter_u: %d\n', iter_u)

    % --- Core velocity solver ---
    u_s = solver_u_core(visc_s, visc, AGlen_s, para);

    % --- Begin Picard iteration ---
    if iter_u > 2
        u_s_now = u_s;

        % preliminary correction vector
        Cs = u_s_now - u_s_lst;

        denom = sqrt(sum(Cs(:).^2) * sum(C(:).^2)); % accurate in math
        denom = max(denom, 1e-20);

        cos_diff_u = dot(Cs(:), C(:)) / denom;
        cos_diff_u = max(min(cos_diff_u, 1.0), -1.0); % numerical safety

        Sita = acos(cos_diff_u); % Sita is within [0, pi]

        % Determine adaptive factor
        if Sita <= pi/8
            mu1 = 2.5;
        elseif Sita >= 19*pi/20
            mu1 = 0.5;
        else
            mu1 = 1.0;
        end

        % relaxed update
        u_s = u_s_lst + mu1 * Cs;

        if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now) < 1e-4
            break
        end     
    end

    if iter_u > 1
        % actual correction vector (both u_s and u_s_lst are relaxed)
        C = u_s - u_s_lst;
    end

    u_s_lst = u_s;

    % --- Update viscosity and strain heating ---
    u = staggerX2main(u_s);
    [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, para);
end
end