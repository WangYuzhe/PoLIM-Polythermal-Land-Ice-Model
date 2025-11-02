function [u_s, u, visc_s, visc, strainHeat] = solver_u_iter_chatgpt(visc_s,visc,AGlen_s,para)
% Wrapped Picard iteration for velocity solver
% Created on 2023/12/1
% Updated on 2025/11/1

global iter_u u_s_lst

iter_u = 0;

while true
    iter_u = iter_u + 1;
    fprintf('iter_u: %d\n', iter_u)

    % --- Core velocity solver ---
    u_s = solver_u_core(visc_s, visc, AGlen_s, para);

    % --- Begin Picard iteration ---
    if iter_u > 2
        u_s_now = u_s;
        Cs = u_s_now - u_s_lst;

        % avoid divide-by-zero
        % ratio_C = dot(Cs(:), C(:)) / (sumsqr(Cs) * sumsqr(C) + 1e-10);

        nume = dot(Cs(:), C(:));
        denom = sqrt(sum(Cs(:).^2) * sum(C(:).^2));
        ratio_C = nume / (denom + 1e-10);

        ratio_C = max(min(ratio_C, 1.0), -1.0); % numerical safety

        Sita = acos(ratio_C);

        % Determine adaptive factor
        if Sita <= pi/8
            mu1 = 2.5;
        elseif Sita >= 19*pi/20
            mu1 = 0.5;
        else
            mu1 = 1.0;
        end

        u_s = u_s_lst + mu1 * Cs;

        % diff_u_now_lst = u_s_now - u_s_lst;
        % ratio_diff_u = sqrt(sum(diff_u_now_lst(:).^2) ./ sum(u_s_now(:).^2));
        % if ratio_diff_u < 1e-4
        %     break
        % end

        if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now) < 1e-4
            break
        end
    end

    if iter_u > 1
        C = u_s - u_s_lst;
    end

    u_s_lst = u_s;

    if iter_u >= para.iter_max
        break
    end

    % --- Update viscosity and strain heating ---
    u = staggerX2main(u_s);
    [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, para);
end
end