function func_main(L)
global SPY N Ms iter iter_max u_s_lst experiment
set_ice_parameters;
[~, hB] = func_set_ice_geometry(L);

% Initialization
AGlen_s = zeros(N, Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,Ms) + 1e13/SPY; % [Pa a]

iter= 0;
while 1
    iter = iter + 1;
    fprintf('iter: %d\n', iter)
    
    u_s = solver_u_period(visc_s, visc);
    
    %---------------------Picard iteration---------------------
    if iter > 2
        u_s_now = u_s;
        Cs = u_s_now - u_s_lst;
        Sita = acos(Cs'*C/(sumsqr(Cs)*sumsqr(C)));
        if isequal(Sita <= pi/8, ones(Ms,Ms))
            mu1 = 2.5;
        elseif isequal(Sita > pi/8,ones(Ms,Ms)) && isequal(Sita < 19*pi/20, ones(Ms,Ms))
            mu1 = 1;
        elseif isequal(Sita >= 19*pi/20, ones(Ms,Ms))
            mu1 = 0.5;
        end
        u_s = u_s_lst + mu1*Cs;
        
        if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now) < 1e-6
            break
        end
    end
    
    if iter > 1
        C = u_s - u_s_lst;
    end
    
    if iter >= iter_max
        break
    end
    %---------------------Picard iteration---------------------
    
    u_s_lst = u_s;
    
    u1 = (u_s(:,1:end-1)+u_s(:,2:end))/2;
    u = [u1(:,end-1), u1];
    
    [visc_s, visc] = get_ice_viscosity(u_s, u, AGlen_s);
end

% calculation
[~, w] = get_ice_w(u_s, u);
[tauxz, ~] = calc_tauxz(u_s, visc_s);
[Pdiff, ~] = calc_Pdiff(u_s, visc_s);

% print the result
fprintf('Max surface velocity: %3.2f \n', max(u(end,1:end-1)))
fprintf('Mean surface velocity: %3.2f \n', mean(u(end,1:end-1)))
%--------------------------------------------------------------------------
% construct the file name

global xi H zeta
xx = ones(N,1)*xi;
yy = ones(N,1)*hB(1:end-2) + zeta*H(1:end-2);
result_u = {xx; yy; u};

result_ismip = [xi/xi(end); u(end,1:end-1); w(end,1:end-1); tauxz; Pdiff];
result_ismip = result_ismip';

if experiment == 1 % Exp. B
    output0 = strcat('wyz1', 'b', sprintf('%03d',L/1000));
    csvwrite(strcat('./results_csv/', output0, '.csv'), result_ismip)
elseif experiment == 2 % Exp. D
    output0 = strcat('wyz1', 'd', sprintf('%03d',L/1000));
    csvwrite(strcat('./results_csv/', output0, '.csv'), result_ismip)
end

end