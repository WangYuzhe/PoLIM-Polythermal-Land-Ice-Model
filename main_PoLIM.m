clc
clearvars
clearvars -global

global SPY N M Ms iter_max iter_u u_s_lst iTimeStep type_thermal_model...
    At_E At_T At_omega At_CTS At_Kappa_s At_Hw At_Ht At_isTemperate...
    At_temperateWaterFluxDarcy

set_ice_parameters;
set_ice_geometry();

dt = 1;
endTime = 50;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

AGlen_s = zeros(N,Ms) + 1e-16;
visc_s = zeros(N,Ms) + 1e13/SPY;
visc = zeros(N,M) + 1e13/SPY;

At_u = zeros(N,M,numTimeStep);
At_w = zeros(N,M,numTimeStep);
At_E                       = zeros(N, M, numTimeStep);
At_T                       = zeros(N, M, numTimeStep);
At_omega                   = zeros(N, M, numTimeStep);
At_Kappa_s                 = zeros(N-1, M, numTimeStep);
At_CTS                     = zeros(numTimeStep, M);
At_Ht                      = zeros(numTimeStep, M);
At_basalMeltRate           = zeros(numTimeStep, M);
At_Hw                      = zeros(numTimeStep, M);
At_basalT                  = zeros(numTimeStep, M);
At_isTemperate             = zeros(numTimeStep, M);
At_temperateWaterFlux      = zeros(N-1, M, numTimeStep);
At_temperateWaterFluxDarcy = zeros(N-1,M,numTimeStep);
At_drainToBed              = zeros(numTimeStep, M);
At_hS = zeros(numTimeStep,M);
At_hB = zeros(numTimeStep,M);
At_H = zeros(numTimeStep,M);

if type_thermal_model ~=3
    Esbc = set_thermalSBC();

    Eini = get_initial_enthalpy(Esbc);
end

is_auto_thermalBasalBC = 1;
type_thermalBasalBC = 1;
has_Greve_drainage = 1;

for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    set_staggered_grid();
    set_thermalSBC();
    
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
        
        [w_vs, w] = get_ice_w(u_s, u);
         
        [~, ~, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s);
        
        switch type_thermal_model
            case 1 % standard enthalpy gradient model
                [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, drainToBed] = ...
                    solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt,...
                    Esbc, Eini, is_auto_thermalBasalBC, type_thermalBasalBC, has_Greve_drainage);
            case 2  % modified enthalpy gradient model
                [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, temperateWaterFluxDarcy] =...
                    solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, dt,...
                    Esbc, Eini, is_auto_thermalBasalBC, type_thermalBasalBC);
            case 3 % isothermal assumption
        end
        
        if type_thermal_model ~= 3
            AGlen_s = get_AGlen(T);
        end

        [visc_s, visc, ~] = get_ice_viscosity(u_s, u, AGlen_s);
    end
    
    At_u(:,:,iTimeStep) = u;
    At_w(:,:,iTimeStep) = w;
    %     At_hS(iTimeStep, :) = hS;
    %     At_hB(iTimeStep, :) = hB;
    %     At_H(iTimeStep, :) = H;
    
    if type_thermal_model ~= 3
        At_E(:,:,iTimeStep) = E;
        At_T(:,:,iTimeStep) = T;
        At_omega(:,:,iTimeStep) = omega;
        At_Kappa_s(:,:,iTimeStep) = Kappa_s;
        At_CTS(iTimeStep,:) = CTS;
        At_Ht(iTimeStep,:) = Ht;
        At_basalMeltRate(iTimeStep,:) = basalMeltRate;
        At_basalT(iTimeStep,:) = T(1,:);
        At_temperateWaterFlux(:,:,iTimeStep) = temperateWaterFlux;

        if iTimeStep == 1
            At_Hw(iTimeStep,:) = zeros(1,M) + dt*basalMeltRate;
        else
            At_Hw(iTimeStep,:) = At_Hw(iTimeStep-1,:) + dt*basalMeltRate;
        end

        logic1 = (At_Hw(iTimeStep,:)>0 & At_isTemperate(iTimeStep,:)==0);
        At_isTemperate(iTimeStep,:) = At_isTemperate(iTimeStep,:) | logic1;
    end
    
    if type_thermal_model==1
        At_drainToBed(iTimeStep,:) = drainToBed;
    end
    
    if type_thermal_model==2
        At_temperateWaterFluxDarcy(:,:,iTimeStep) = temperateWaterFluxDarcy;
    end
    
    
    fprintf('Mean surface velocity: %3.2f \n', mean(u(end,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
end

plot_enthalpy_field
% plot_uField_uSurf
