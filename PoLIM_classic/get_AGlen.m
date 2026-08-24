function [AGlen_s] = get_AGlen(T, omega, CTS, p)
% Date: 2017-12-17
% Author: Wang Yuzhe
% Calculate the flow rate factor A.

% Input:
% T: temperature
% para: parameters (struct)

% Output:
% AGlen_s: flow rate factor A on staggered grid.

global N M H zeta

SPY = p.SPY;
rhoi = p.rhoi;
g = p.g;
R = p.R;
betaCC = p.betaCC;
type_Arrhenius = p.type_Arrhenius;
is_duval = p.is_duval;

AGlen = zeros(N,M);
if strcmpi(type_Arrhenius,'greve')
    AGlen0_cold = 3.985e-13;
    AGlen0_warm = 1.916e3;
    
    Q_cold = 60e3;
    Q_warm = 139e3;
    for i = 1:M
        T_Corr = T(:,i) + betaCC*rhoi*g*(1-zeta).*H(i);
        idx_le_n10 = (T_Corr <= 263.15);
        idx_ge_n10 = (T_Corr > 263.15);
        
        AGlen(:,i) = AGlen0_cold*exp(-Q_cold./(R*T_Corr)).*idx_le_n10 +...
            AGlen0_warm*exp(-Q_warm./(R*T_Corr)).*idx_ge_n10; % [Pa-3 s-1]
    end
    
elseif strcmpi(type_Arrhenius,'cuffey')
    AGlen0 = 3.5e-25; % constant prefactor at -10 celsius [Pa-3 s-1], Cuffey&Paterson (2010), p74, eq3.36
    Q_cold = 60e3; % activation energy for creep [J mol-1]
    Q_warm = 115e3; % activation energy for creep [J mol-1]

    % loop from the upstream to the downstream columns
    for i = 1:M
        delta_Tpmp = betaCC*rhoi*g*H(i)*(1-zeta);
        Tref_corr = 263.15 + delta_Tpmp;
        T_corr = T(:,i) + delta_Tpmp;
        iCTS = CTS(i);

        % cold layer
        idx_COLD = (iCTS+1):N;
        idx_T_le_m10C = (T(idx_COLD,i)<=263.15);
        idx_T_gt_m10C = (T(idx_COLD,i)>263.15);
        T_corr_COLD = T_corr(idx_COLD);
        Tref_corr_COLD = Tref_corr(idx_COLD);
        AGlen(idx_COLD,i) = AGlen0*exp(-Q_cold./R.*(1./T_corr_COLD-1./Tref_corr_COLD)).*idx_T_le_m10C +...
            AGlen0*exp(-Q_warm./R.*(1./T_corr_COLD-1./Tref_corr_COLD)).*idx_T_gt_m10C; % [Pa-3 s-1]

        % temperate layer
        if iCTS>=1
            % index for temperate ice
            idx_TEMP = 1:iCTS;

            % set the water content below the upper bound (3%)
            idx_le_omega_limit = (omega(:,i)<=p.omega_max);
            idx_gt_omega_limit = (omega(:,i)>p.omega_max);
            omega(:,i) = omega(:,i).*idx_le_omega_limit + p.omega_max*idx_gt_omega_limit;

            % Duval-Lliboutry creep relation for temperate ice
            AGlen(idx_TEMP,i) = 24e-25*(1 + p.is_duval*181.25*omega(idx_TEMP,i)); % [Pa-1 s-1]
        end
    end
else
    disp('type_Arrhenius should be "greve", "cuffey"')
end

AGlen = AGlen*SPY; % [Pa-3 a-1]

AGlen_s = main2staggerX(AGlen);
end