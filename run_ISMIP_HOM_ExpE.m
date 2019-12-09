% Run the ISMIP-HOM experiment E (Exp. E: Haut Glacier d¡¯Arolla)
% id_exp = 'e000': the first experiment (E1)
% id_exp = 'e001': the second experiment (E2)

clc
clearvars
clearvars -global

global SPY N M Ms iter_max iter_u u_s_lst iTimeStep

id_exp = 'e000'; % expE2: e001
if strcmp(id_exp, 'e000')
    set_ice_parameters_ismip_hom(1,2,3,0) % id_type_BBC, id_type_geometry, id_type_valley, id_isFlowband
else
    set_ice_parameters_ismip_hom(4,2,3,0)
end

set_ice_geometry();

% Initialization
AGlen_s = zeros(N, Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N, Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N, M) + 1e13/SPY; % [Pa a]

numTimeStep = 1;
At_u = zeros(N, M, numTimeStep);
At_w = zeros(N, M, numTimeStep);
At_hS = zeros(numTimeStep, M);
At_hB = zeros(numTimeStep, M);
At_H = zeros(numTimeStep, M);

for iTimeStep = 1:numTimeStep
    
    fprintf('iTimeStep: %d \n', iTimeStep)
    set_staggered_grid();
    
    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        fprintf('iter_u: %d\n', iter_u)
        
        [u_s] = solver_u(visc_s, visc, AGlen_s);
        
        %-------------------------Picard iteration-------------------------
        if iter_u > 2
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
            
            if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now) < 1e-4
                break
            end
        end
        
        if iter_u > 1
            C = u_s - u_s_lst;
        end
        
        if iter_u >= iter_max
            break
        end
        u_s_lst = u_s;
        %-------------------------Picard iteration-------------------------
        
        u = staggerX2main(u_s);
        
        % update the viscosity
        [visc_s, visc, ~] = get_ice_viscosity(u_s, u, AGlen_s);
    end
    
    At_u(:,:,iTimeStep) = u;
    
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
end

[w_vs, w] = get_ice_w(u_s, u);

[tauxz, tauxz_s] = calc_tauxz(u_s, visc_s);
[Pdiff, Pdiff_s] = calc_Pdiff(u_s, visc_s);

% Figure
plot_uField_uSurf

% Save results
if strcmp(id_exp, 'e000')
    filename = strcat('wyz1', 'e000');
    filename1 = strcat('fullwyz1', 'e000');
    
    result = [xi/xi(end); u(end,:); w(end,:); tauxz; Pdiff];
    result = result';
    
    result1 = {xi/xi(end); u(end,:); w(end,:); tauxz; Pdiff; xx; yy; u};
    
    eval([filename '=' 'result' ';'])
    eval('save(filename, filename)')
    
    eval([filename1 '=' 'result1' ';'])
    eval('save(filename1, filename1)')
    
else
    filename = strcat('wyz1', 'e001');
    filename1 = strcat('fullwyz1', 'e001');
    
    result = [xi/xi(end); u(end,:); w(end,:); tauxz; Pdiff];
    result = result';
    
    result1 = {xi/xi(end); u(end,:); w(end,:); tauxz; Pdiff; xx; yy; u};
    
    eval([filename '=' 'result' ';'])
    eval('save(filename, filename)')
    
    eval([filename1 '=' 'result1' ';'])
    eval('save(filename1, filename1)')
end