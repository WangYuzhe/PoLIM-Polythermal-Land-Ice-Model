% Date: 2017-7-12
% Author: Wang Yuzhe
% Calculate the glacier thermal regime using the enthalpy scheme.
% Use the enthalpy definition and the water drainage model proposed by Anderas Aschwanden.

function [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, drainToBed] = ...
    solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc, Eini, is_auto_thermalBasalBC, type_thermalBasalBC, has_Greve_drainage)
global SPY rho g kc Cp Qgeo betaCC rhow Lw Tref Kc Kt...
    M N dx dzeta zeta H dzetadx iTimeStep...
    At_E At_Hw At_Ht At_Kappa_s At_isTemperate At_omega

dt = dt*SPY;
u = u/SPY;
u_s = u_s/SPY;
w = w/SPY;
w_vs = w_vs/SPY;
strainHeat = strainHeat/SPY;

dzetadx_vs = (dzetadx(1:end-1,:) + dzetadx(2:end,:))/2;
u_vs = (u(1:end-1,:) + u(2:end,:))/2;

coeff = u.*dzetadx + w./(ones(N,1)*H);
coeff_vs = u_vs.*dzetadx_vs + w_vs(1:N-1, :)./(ones(N-1,1)*H);

LT1 = zeros(N,1);
LT2 = zeros(N,1);
LT3 = zeros(N,1);
RT = zeros(N,1);

E = zeros(N,M);
T = zeros(N,M);
omega = zeros(N,M);

Ht = zeros(1,M);
basalMeltRate = zeros(1,M);
Kappa = zeros(N,1);
CTS = zeros(1,M);
temperateWaterFlux = zeros(N-1,M);
drainageRate = zeros(N,1);
drainToBed = zeros(1,M);

if iTimeStep == 1
    isTransient = 0;
    Enm1 = Eini;
    Kappa_s = Kc*ones(N-1,M);
    omega4drainage = zeros(N,M);
else
    isTransient = 1;
    Enm1 = At_E(:,:,iTimeStep-1);
    Kappa_s = At_Kappa_s(:,:,iTimeStep-1);
    omega4drainage = At_omega(:,:,iTimeStep-1);
end

for i = 1:M
        
    Tpmp_i = 273.15 - betaCC*rho*g*H(i)*(1-zeta);
    Epmp_i = Cp*(Tpmp_i - Tref);
    
    for j = 2:N-1
        if has_Greve_drainage
            drainageRate(j) = drainageFunc(omega4drainage(j,i))/SPY;
        else
            drainageRate(j) = 0;
        end
        
        if i == 1
            RT(j) = isTransient*Enm1(j,i)/(dt) + strainHeat(j,i)/rho;
        else
            RT(j) = isTransient*Enm1(j,i)/(dt) + strainHeat(j,i)/rho...
                + u_s(j,i)*Enm1(j,i-1)/dx...
				- rhow/rho*Lw*drainageRate(j);
        end
        
        if coeff(j,i) > 0
            LT1(j) = -coeff_vs(j-1,i)/dzeta - Kappa_s(j-1,i)/(H(i)^2*dzeta^2);
            LT2(j) = isTransient/(dt) + u_s(j,i)/dx + coeff_vs(j-1,i)/dzeta + (Kappa_s(j-1,i)+Kappa_s(j,i))/(H(i)^2*dzeta^2);
            LT3(j) = -Kappa_s(j,i)/(H(i)^2*dzeta^2);
        elseif coeff(j,i) <= 0
            LT1(j) = -Kappa_s(j-1,i)/(H(i)^2*dzeta^2);
            LT2(j) = isTransient/(dt) + u_s(j,i)/dx - coeff_vs(j,i)/dzeta + (Kappa_s(j-1,i)+Kappa_s(j,i))/(H(i)^2*dzeta^2);
            LT3(j) = coeff_vs(j,i)/dzeta - Kappa_s(j,i)/(H(i)^2*dzeta^2);
        end
    end

    LT1(N) = 0; LT2(N) = 1; LT3(N) = 0; RT(N) = Esbc(i);

    if is_auto_thermalBasalBC
        if iTimeStep == 1
            [LT, RT] = basalBC_cold_base_dry(LT1, LT2, LT3, RT, H(i));
        else
            if At_isTemperate(iTimeStep-1,i) == 0
                if At_Hw(iTimeStep-1,i) > 0
                    [LT, RT] = basalBC_cold_base_wet(LT1, LT2, LT3, RT, Epmp_i);
                else
                    [LT, RT] = basalBC_cold_base_dry(LT1, LT2, LT3, RT, H(i));
                end
            else
                if At_Ht(iTimeStep-1,i) > 0
                    [LT, RT] = basalBC_temperate_layer(LT1, LT2, LT3, RT);
                else
                    [LT, RT] = basalBC_temperate_base(LT1, LT2, LT3, RT, Epmp_i);
                end
            end
        end
    else
        switch type_thermalBasalBC
            case 1
                [LT, RT] = basalBC_cold_base_dry(LT1, LT2, LT3, RT, H(i));
            case 2
                [LT, RT] = basalBC_temperate_base(LT1, LT2, LT3, RT, Epmp_i);
            case 3
                [LT, RT] = basalBC_temperate_layer(LT1, LT2, LT3, RT);
            case 4
                [LT, RT] = basalBC_cold_base_wet(LT1, LT2, LT3, RT, Epmp_i);
        end        
    end

    E(:,i) = LT\RT;
    
    logic1 = E(:,i) >= Epmp_i;
    T(:,i) = logic1.*Tpmp_i + (1-logic1).*(E(:,i)/Cp + Tref);
    omega(:,i) = logic1.*(E(:,i) - Epmp_i)/Lw;

    jcts = find(logic1, 1, 'last');
    
    if isempty(jcts)
        CTS(i) = 0;
        Ht(i) = 0;
        basalMeltRate(i) = 0;
    else
        CTS(i) = jcts;
        Ht(i) = (jcts-1)*H(i)*dzeta;
        Kappa(1:jcts-1) = Kt;
        Kappa(jcts:end) = Kc;
        Kappa_harmmean = harmmean([Kappa(1:end-1)'; Kappa(2:end)']);
        Kappa_s(:,i) = Kappa_harmmean';
        
        if is_auto_thermalBasalBC
            if jcts == 1
                [LT, RT] = basalBC_temperate_base(LT1, LT2, LT3, RT, Epmp_i);
            else
                [LT, RT] = basalBC_temperate_layer(LT1, LT2, LT3, RT);
            end
        else
            switch type_thermalBasalBC
                case 1
                    [LT, RT] = basalBC_cold_base_dry(LT1, LT2, LT3, RT, H(i));
                case 2
                    [LT, RT] = basalBC_temperate_base(LT1, LT2, LT3, RT, Epmp_i);
                case 3
                    [LT, RT] = basalBC_temperate_layer(LT1, LT2, LT3, RT);
                case 4
                    [LT, RT] = basalBC_cold_base_wet(LT1, LT2, LT3, RT, Epmp_i);
            end
        end
        
        E(:,i) = LT\RT;
        basalMeltRate(i) = (Qgeo + kc/Cp*(E(2,i)-E(1,i))/(H(i)*dzeta))/(rhow*Lw)*SPY; % [m a-1]
    end

    temperateWaterFlux(:,i) = -Kt*(omega(2:end,i) - omega(1:end-1,i))/(H(i)*dzeta);

    drainToBed(i) = sum(drainageRate)*H(i)*dzeta;

    index1 = omega(:,i) > 3/100;
    index2 = omega(:,i) <= 3/100;
    omega(:,i) = index1*3/100 + index2.*omega(:,i);
end

end
