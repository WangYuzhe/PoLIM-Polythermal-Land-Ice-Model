% Date: 2017-7-12
% Author: Wang Yuzhe
% Calculate the glacier thermal regime using the enthalpy scheme.
% Use the enthalpy definition by Andreas Aschwanden and use the water drainage model proposed by Ian Hewitt

function [E, T, omega, Kappa_s, CTS, Ht, basalMeltRate, temperateWaterFlux, temperateWaterFluxDarcy] =...
    solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, dt, Esbc, Eini, is_auto_thermalBasalBC, type_thermalBasalBC)
global SPY rho g kc Cp Qgeo betaCC rhow Lw Tref Kc Kt k0 eta_w
global M N dx dzeta zeta H dzetadx
global iTimeStep At_E At_Hw At_Ht At_Kappa_s At_isTemperate At_temperateWaterFluxDarcy

% convert [xx a-1] to [xx s-1]
dt = dt*SPY;
u = u/SPY; % [m s-1]
u_s = u_s/SPY; % [m s-1]
w = w/SPY;  % [m s-1]
w_vs = w_vs/SPY; % [m s-1]
strainHeat = strainHeat/SPY; % [Pa s-1];

% 'var_vs' means the secondary grid point in zeta coordinate
dzetadx_vs = (dzetadx(1:end-1,:) + dzetadx(2:end,:))/2; % <N-1 * M>
u_vs = (u(1:end-1,:) + u(2:end,:))/2; % <N-1 * M>

% coefficient for the vertical derivative of enthalpy (\partial E / \partial zeta)
coeff = u.*dzetadx + w./(ones(N,1)*H); % [s-1]; <N * M>
coeff_vs = u_vs.*dzetadx_vs + w_vs(1:N-1, :)./(ones(N-1,1)*H); % [s-1]; <N-1 * M>

% initialization
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
temperateWaterFluxDarcy = zeros(N-1,M);
temperateWaterFluxDiffusive = zeros(N-1,M);

if iTimeStep == 1
    isTransient = 0;
    Enm1 = Eini;
    Kappa_s = Kc*ones(N-1,M);    
    temperateWaterFluxDarcy_last = zeros(N,M);
else
    isTransient = 1;
    Enm1 = At_E(:,:,iTimeStep-1);
    Kappa_s = At_Kappa_s(:,:,iTimeStep-1);
    temperateWaterFluxDarcy_last = At_temperateWaterFluxDarcy(:,:,iTimeStep-1);
end

for i = 1:M
    
    % Enthalpy at the pressure-melting point for the i-th column
    Tpmp_i = 273.15 - betaCC*rho*g*H(i)*(1-zeta); % <N * 1>
    Epmp_i = Cp*(Tpmp_i - Tref); % <N * 1>
    
    for j = 2:N-1
        if i == 1
            RT(j) = isTransient*Enm1(j,i)/(dt) + strainHeat(j,i)/rho;
        else
            RT(j) = isTransient*Enm1(j,i)/(dt) + strainHeat(j,i)/rho...
                + u_s(j,i)*Enm1(j,i-1)/dx...
				- rhow/rho*Lw*(temperateWaterFluxDarcy_last(j,i)-temperateWaterFluxDarcy_last(j-1,i))/(H(i)*dzeta);
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
    
    % surface BC
    LT1(N) = 0; LT2(N) = 1; LT3(N) = 0; RT(N) = Esbc(i);
    
    % basal BC
    if is_auto_thermalBasalBC % decision chart (Aschwanden et al., 2012, Figure 5)
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
    else % Do not use the decision chart. But specify a type of basal boundary condition.
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
    
    % solution
    E(:,i) = LT\RT; % <N * 1>
    
    logic1 = E(:,i) >= Epmp_i; % 1: temperate; 0: cold. <N * 1>
    T(:,i) = logic1.*Tpmp_i + (1-logic1).*(E(:,i)/Cp + Tref);
    omega(:,i) = logic1.*(E(:,i) - Epmp_i)/Lw;
        
    % If cold base, 'jcts' is an empty array.
    % If temperate base or temperate layer, 'jcts' >= 1
    jcts = find(logic1, 1, 'last');
    
    if isempty(jcts)      %% COLD
        CTS(i) = 0;
        Ht(i) = 0;
        basalMeltRate(i) = 0;
    else                  %% TEMPERATE
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
    
    % calculate the water flux
%     omega_i_vs = (omega(2:end,i)+omega(1:end-1,i))/2;
%     temperateWaterFlux(:,i) = -k0*omega_i_vs.^2*(rhow-rho)*g/(eta_w);
%     temperateWaterFlux(:,i) = -k0*omega(:,i).^2*(rhow-rho)*g/(eta_w);
%     temperateWaterFlux(:,i) = -k0*omega(2:end,i).^2*(rhow-rho)*g/(eta_w);
    temperateWaterFluxDiffusive(:,i) = -Kt*(omega(2:end,i) - omega(1:end-1,i))/(H(i)*dzeta);
    temperateWaterFluxDarcy(:,i) = -k0*omega(2:end,i).^2*(rhow-rho)*g/(eta_w);
    temperateWaterFlux(:,i) = temperateWaterFluxDiffusive(:,i) + temperateWaterFluxDarcy(:,i);

    % set the upper bound of the water content (3%)
    index1 = omega(:,i) > 3/100;
    index2 = omega(:,i) <= 3/100;
    omega(:,i) = index1*3/100 + index2.*omega(:,i);
end

end