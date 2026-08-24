function [E, T, omega, Kappa_vs, CTS, is_TEMP, thk_TEMP, thk_w, m_basal, qw_TEMP_drain_greve, qw_TEMP_diffu] = ...
    solve_enth_SEGM(u, u_s, w, w_vs, strain_heat, friction_heat, dt, Esbc, Eini, p)
% Date: 2017-7-12
% Author: Wang Yuzhe
% Calculate the glacier enthalpy field
% Use the enthalpy definition by Aschwanden et al. (2012)
% Use the water drainage model proposed by Greve (2009)
% Use the one-layer melting CTS scheme proposed by Blatter&Greve (2015)
% <is_TEMP>, <thk_TEMP>, and <thk_w> determine the basal boundary condition

global M N dx dzeta zeta H dzetadx iTimeStep enth_lst

%% decision chart for basal location to determine boundary condition
    function [LT, RT] = basalBC_COLD_base_dry(LT1, LT2, LT3, RT, H_icol, para)
        LT1(1) = 0;
        LT2(1) = -1/dzeta;
        LT3(1) = 1/dzeta;
        RT(1) = -para.Qgeo*H_icol*para.Cp/para.kc;
        LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
    end

    function [LT, RT] = basalBC_COLD_base_wet(LT1, LT2, LT3, RT, Epmp_i)
        LT1(1) = 0;
        LT2(1) = 1;
        LT3(1) = 0;
        RT(1) = Epmp_i(1);
        LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
    end

    function [LT, RT] = basalBC_TEMP_base(LT1, LT2, LT3, RT, Epmp_i)
        LT1(1) = 0;
        LT2(1) = 1;
        LT3(1) = 0;
        RT(1) = Epmp_i(1);

        % RT(1) = Epmp_i(1) + 3.34e5*0.01;

        LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
    end

    function [LT, RT] = basalBC_TEMP_layer(LT1, LT2, LT3, RT)
        LT1(1) = 0;
        LT2(1) = -2;
        LT3(1) = 3;
        RT(1) = 0;
        LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
        LT(1,3) = -1;

        % LT1(1) = 0;
        % LT2(1) = -1;
        % LT3(1) = 1;
        % RT(1) = 0;
        % LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
    end

%% parameters
SPY = p.SPY;
rhoi = p.rhoi;
rhow = p.rhow;
g = p.g;
kc = p.kc;
Cp = p.Cp;
Lw = p.Lw;
Tref = p.Tref;
Kc = p.Kc;
Kt = p.Kt;
is_enth_trans = p.is_enth_trans; % transient solver

%% initialization
% convert [xx a-1] to [xx s-1]
dt = dt*SPY; % [s-1]
u = u/SPY; % [m s-1]
u_s = u_s/SPY; % [m s-1]
w = w/SPY;  % [m s-1]
w_vs = w_vs/SPY; % [m s-1]
strain_heat = strain_heat/SPY; % [Pa s-1];

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

is_TEMP = zeros(1,M);
thk_TEMP = zeros(1,M);
thk_w = zeros(1,M);
m_basal = zeros(1,M);
CTS = zeros(1,M);
Kappa_i = zeros(N,1);
Kappa_vs = zeros(N-1,M);
qw_TEMP_diffu = zeros(N-1,M);
qw_TEMP_drain_greve = zeros(1,M);
rate_TEMP_water = zeros(N,1);

if iTimeStep==1
    Elst = Eini;
    Kappa_vs_lst = Kc*ones(N-1,M);
    omega_lst = zeros(N,M);
    is_TEMP_lst = zeros(1,M);
    thk_TEMP_lst = zeros(1,M);
    thk_w_lst = zeros(1,M);
else
    Elst = enth_lst.E;
    Kappa_vs_lst = enth_lst.Kappa_vs;
    omega_lst = enth_lst.omega;
    is_TEMP_lst = enth_lst.is_TEMP;
    thk_TEMP_lst = enth_lst.thk_TEMP;
    thk_w_lst = enth_lst.thk_w;
end

for i=1:M
    % Enthalpy at the pressure-melting point for the i-th column
    H_i = H(i);
    Tpmp_i = 273.15 - p.betaCC*rhoi*g*H_i*(1-zeta); % <N * 1>
    Epmp_i = Cp*(Tpmp_i - Tref); % <N * 1>
    
    for j = 2:N-1
        if p.is_greve_drain
            rate_TEMP_water(j) = drain_temp_water_greve(omega_lst(j,i))/SPY; % [s-1]
        else
            rate_TEMP_water(j) = 0.0;
        end
        
        if M==1 % ice slab (no horizontal advection and diffusion)
            RT(j) = is_enth_trans*Elst(j,i)/dt + strain_heat(j,i)/rhoi -...
                rhow/rhoi*Lw*rate_TEMP_water(j);
        else
            if i==1
                RT(j) = is_enth_trans*Elst(j,i)/dt + strain_heat(j,i)/rhoi;
            else
                RT(j) = is_enth_trans*Elst(j,i)/dt + strain_heat(j,i)/rhoi -...
                    rhow/rhoi*Lw*rate_TEMP_water(j) +...
                    u_s(j,i)*Elst(j,i-1)/dx;
            end
        end
        
        if coeff(j,i)>0
            LT1(j) = -coeff_vs(j-1,i)/dzeta - Kappa_vs_lst(j-1,i)/(H_i^2*dzeta^2);
            LT2(j) = is_enth_trans/dt + u_s(j,i)/dx + coeff_vs(j-1,i)/dzeta +...
                (Kappa_vs_lst(j-1,i)+Kappa_vs_lst(j,i))/(H_i^2*dzeta^2);
            LT3(j) = -Kappa_vs_lst(j,i)/(H_i^2*dzeta^2);
        elseif coeff(j,i)<=0
            LT1(j) = -Kappa_vs_lst(j-1,i)/(H_i^2*dzeta^2);
            LT2(j) = is_enth_trans/dt + u_s(j,i)/dx - coeff_vs(j,i)/dzeta +...
                (Kappa_vs_lst(j-1,i)+Kappa_vs_lst(j,i))/(H_i^2*dzeta^2);
            LT3(j) = coeff_vs(j,i)/dzeta - Kappa_vs_lst(j,i)/(H_i^2*dzeta^2);
        end
    end

    % surface BC
    LT1(N) = 0;
    LT2(N) = 1;
    LT3(N) = 0;
    RT(N) = Esbc(i);

    % basal BC
    if p.is_auto_enth_BBC % decision chart (Aschwanden et al., 2012, Figure 5)
        if iTimeStep==1
            [LT, RT] = basalBC_COLD_base_dry(LT1, LT2, LT3, RT, H_i, p);
        else
            if is_TEMP_lst(i)==0
                if thk_w_lst(i)>0
                    [LT, RT] = basalBC_COLD_base_wet(LT1, LT2, LT3, RT, Epmp_i);
                else
                    [LT, RT] = basalBC_COLD_base_dry(LT1, LT2, LT3, RT, H_i, p);
                end
            else
                if thk_TEMP_lst(i)>0
                    [LT, RT] = basalBC_TEMP_layer(LT1, LT2, LT3, RT);
                else
                    [LT, RT] = basalBC_TEMP_base(LT1, LT2, LT3, RT, Epmp_i);
                end
            end
        end
    else % Specify a type of basal boundary condition
        switch p.type_enth_BBC
            case 'COLD_base_dry'
                [LT, RT] = basalBC_COLD_base_dry(LT1, LT2, LT3, RT, H_i, p);
            case 'TEMP_base'
                [LT, RT] = basalBC_TEMP_base(LT1, LT2, LT3, RT, Epmp_i);
            case 'TEMP_layer'
                [LT, RT] = basalBC_TEMP_layer(LT1, LT2, LT3, RT);
            case 'COLD_base_wet'
                [LT, RT] = basalBC_COLD_base_wet(LT1, LT2, LT3, RT, Epmp_i);
        end
    end

    % solution of the predictor step
    E(:,i) = LT\RT; % <N * 1>

    logic1 = E(:,i)>=Epmp_i; % 1: temperate; 0: cold. <N * 1>
    jcts = find(logic1, 1, 'last'); % empty: cold; jcts>=1: temperate

    if isempty(jcts) %% COLD
        CTS(i) = 0;
        thk_TEMP(i) = 0;
        Kappa_vs(:,i) = Kc;

        T(:,i) = E(:,i)/Cp + Tref;
        omega(:,i) = 0;
        m_basal(i) = 0;

    else %% TEMP
        CTS(i) = jcts;
        thk_TEMP(i) = (jcts-1)*H_i*dzeta;

        Kappa_i(1:jcts-1) = Kt;
        Kappa_i(jcts:end) = Kc;
        Kappa_harmmean = harmmean([Kappa_i(1:end-1)'; Kappa_i(2:end)']);
        Kappa_vs(:,i) = Kappa_harmmean'; % harmonic mean

        if p.is_blatter_meltCTS
            % If a temperate layer existed, correct the enthalpy of the cold layer
            % <BEGIN corrector step> only for cold layer (Blatter & Greve, 2015, p201)
            if jcts>1
                Nc = N - jcts + 1;
                LT1c = zeros(Nc, 1);
                LT2c = zeros(Nc, 1);
                LT3c = zeros(Nc, 1);
                RTc = zeros(Nc, 1);
                for jc = 2:(Nc-1)
                    jj = jc + jcts - 1;
                    if M==1 % ice slab case (no horizontal advection and diffusion)
                        RTc(jc) = E(jj,i)/dt + strain_heat(jj,i)/rhoi;
                    else
                        if i==1
                            RTc(jc) = E(jj,i)/dt + strain_heat(jj,i)/rhoi;
                        else
                            RTc(jc) = E(jj,i)/dt + strain_heat(jj,i)/rhoi +...
                                u_s(jj,i)*E(jj,i-1)/dx;
                        end
                    end

                    if coeff(jj,i)>0
                        LT1c(jc) = -coeff_vs(jj-1,i)/dzeta - Kappa_vs(jj-1,i)/(H_i^2*dzeta^2);
                        LT2c(jc) = 1/dt + u_s(jj,i)/dx + coeff_vs(jj-1,i)/dzeta +...
                            (Kappa_vs(jj-1,i)+Kappa_vs(jj,i))/(H_i^2*dzeta^2);
                        LT3c(jc) = -Kappa_vs(jj,i)/(H_i^2*dzeta^2);
                    elseif coeff(jj,i)<=0
                        LT1c(jc) = -Kappa_vs(jj-1,i)/(H_i^2*dzeta^2);
                        LT2c(jc) = 1/dt + u_s(jj,i)/dx - coeff_vs(jj,i)/dzeta +...
                            (Kappa_vs(jj-1,i)+Kappa_vs(jj,i))/(H_i^2*dzeta^2);
                        LT3c(jc) = coeff_vs(jj,i)/dzeta - Kappa_vs(jj,i)/(H_i^2*dzeta^2);
                    end
                end

                % surface BC for cold layer <corrector step>
                LT1c(end) = 0;
                LT2c(end) = 1;
                LT3c(end) = 0;
                RTc(end) = Esbc(i);

                % basal BC for cold layer <corrector step>
                LT1c(1) = 0;
                LT2c(1) = 1;
                LT3c(1) = 0;
                RTc(1) = Epmp_i(jcts);

                LTc = spdiags([[LT1c(2:end);0], LT2c, [0;LT3c(1:end-1)]], [-1,0,1], Nc, Nc);
                Ec = LTc\RTc;
                E(jcts:end,i) = Ec;
            end
            % <END corrector step>
        end

        % update
        T(:,i) = logic1.*Tpmp_i + (1 - logic1).*(E(:,i)/Cp + Tref);
        omega(:,i) = logic1.*(E(:,i) - Epmp_i)/Lw;
        conudct_heat_i = kc*(T(2,i)-T(1,i))/(H_i*dzeta);
        m_basal(i) = (p.Qgeo + conudct_heat_i + friction_heat(i))/(rhow*Lw)*SPY; % basal melt [m a-1]
    end

    thk_w(i) = thk_w_lst(i) + dt/SPY*m_basal(i);

    is_TEMP(i) = (E(1,i)>=Epmp_i(1)) & (thk_w(i)>0);
    
    % diffusive water flux through temperate ice
    qw_TEMP_diffu(:,i) = -Kt*(omega(2:end,i)-omega(1:end-1,i))/(H_i*dzeta);
    
    % vertically integrated drainage [m s-1]
    % drain instantaneously to the bed using the Ralf Greve's function
    qw_TEMP_drain_greve(i) = sum(rate_TEMP_water)*H_i*dzeta; % [m s-1]
    
    % set the upper bound of the water content (3%)
    %     index1 = omega_i > 3/100;
    %     index2 = omega_i <= 3/100;
    %     omega_i = index1*3/100 + index2.*omega_i;
end

end