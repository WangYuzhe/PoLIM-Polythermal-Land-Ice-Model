% Date: 2019-6-17

clear, clc


global SPD SPY type_glacier_geometry M rhow g hB dx

tune_run = 'A3';
set_ice_parameters;

type_glacier_geometry = 2;
% 1: haut d'arolla
% 2: ice sheet
% 3: valley glacier
set_ice_geometry();

% TIME SETTING
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);


% WATER SOURCES
if strcmp(tune_run, 'A3')
    basalMelt = 5.79e-9*ones(1,M);
else
    basalMelt = 4.5e-8*ones(1,M);
end

mouInput = zeros(1,M);
seasonRunoff = zeros(numTimeStep,M);

% PARAMETERS
ub = 1e-6*ones(1,M);
ev = 0;

Ks_range = [1e-6, 1e-5, 1e-4, 2e-4, 4e-4, 5e-4, 1e-3, 1e-2, 1e-1];
visc_ice_b_range = [1e12, 2e12, 3e12, 4e12, 1e13, 5e13, 1e14, 5e14, 1e15];

n_Ks = length(Ks_range);
n_visc_ice_b = length(visc_ice_b_range);

tune_result_N = cell(n_Ks, n_visc_ice_b);
for i = 1:n_Ks
    Ks = Ks_range(i);
    for j = 1:n_visc_ice_b
        fprintf('i = %d, j = %d\n', i, j)
        visc_ice_b = visc_ice_b_range(j);
        
        iter = 0;
        while 1
            iter = iter + 1;
            % INITIAL CONDITIONS
            if iter == 1
                phi_nm1 = rhow*g*hB;
                Hw_nm1 = 1e-1*ones(1,M); % [m]
                dHwdx_nm1 = [(Hw_nm1(2)-Hw_nm1(1))/dx,...
                    (Hw_nm1(3:end) - Hw_nm1(1:end-2))/(2*dx),...
                    (Hw_nm1(end)-Hw_nm1(end-1))/dx];
                dHw3dx_nm1 = 3*Hw_nm1.^2.*dHwdx_nm1;
            end;            
            % SUBGLACIAL HYDROLOGY MODEL
            [Neff, phi, qw, Hw, dHw3dx] = subHydro_storage_Hoffman_tune(dt, basalMelt, mouInput, seasonRunoff, Hw_nm1, dHw3dx_nm1, phi_nm1, ev, ub, Ks, visc_ice_b);
            tune_result_N{i,j} = Neff;
            phi_nm1 = phi;
            Hw_nm1 = Hw;
            dHw3dx_nm1 = dHw3dx;
            if iter > 365
                break
            end
        end
    end
    
end
matname = strcat('tune_result_N_', tune_run);

eval(['save ' matname ' tune_result_N'])