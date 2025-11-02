clc
clearvars
clearvars -global

tic
global N M H dzetadx zeta dzeta iTimeStep...
    At_E At_T At_omega At_CTS At_Kappa_vs...
    At_isTEMP At_thk_TEMP At_thk_w At_qw_TEMP_darcy At_qw_TEMP
%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params_hewitt_slab();
SPY = p.SPY;

p.type_enth_solver = 'MEGM';
p.type_enth_BBC = 2;
p.has_Greve_drainage = 0;
p.k0 = 1e-12;

p.layers = 201;
%% GLACIER GEOMETRY
%
%
set_ice_geometry('./inputs/geo_polythermal_slab', p);
dzetadx = zeros(N,M);
%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 2000; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, endTime);

%% INITIALIZATION
%
%
alpha = 4*pi/180; % inclination angle
vz = -0.2; % [m a-1] vertical velocity

u = zeros(N,M);
u_s = zeros(N,M+1);
w = vz*ones(N,M);
w_vs = vz*ones(N,M);
strainHeat = 2*p.AGlen*(p.rhoi*p.g*sin(alpha)).^4.*H.^4.*(1-zeta).^4*SPY;

%  Related to enthalpy solver
At_E = zeros(N,M,numTimeStep);
At_T = zeros(N,M,numTimeStep);
At_omega = zeros(N,M,numTimeStep);
At_Kappa_vs = zeros(N-1,M,numTimeStep);
At_CTS = zeros(numTimeStep,M);
At_isTEMP = zeros(numTimeStep,M);
At_thk_TEMP = zeros(numTimeStep,M);
At_thk_w = zeros(numTimeStep,M); % water thickness
At_qw_TEMP_darcy = zeros(N-1,M,numTimeStep);
At_qw_TEMP = zeros(N-1,M,numTimeStep);
At_m_basal = zeros(numTimeStep,M); % basal melt rate

% Hewitt-schoof polythermal slab
Esbc = p.Cp*(-1 + 273.15 - p.Tref);
Eini = p.Cp*(-2 + 273.15 - p.Tref)*ones(N,M);

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    % calculate the enthalpy
    if strcmpi(p.type_enth_solver,'SEGM')
        [E, T, omega, Kappa_vs, CTS, thk_TEMP, qw_greveDrain, qw_TEMP_diffu] = ...
            solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt_u, Esbc, Eini, p);
    elseif strcmpi(p.type_enth_solver,'MEGM')
        [E, T, omega, Kappa_vs, CTS, thk_TEMP, qw_TEMP, qw_TEMP_darcy] =...
            solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, dt_u, Esbc, Eini, p);
    elseif strcmpi(p.type_enth_solver,'isoT')
        % pass
    end

    tauxz = zeros(1,M);
    m_basal = (p.Qgeo + p.kc*(T(2,:)-T(1,:))./(H*dzeta) +...
        tauxz.*u(1,:)/SPY)/(p.rhow*p.Lw)*SPY; % basal melt [m a-1]
    m_basal(m_basal<0)=0;
        
    % STORAGE RESULTS   
    At_E(:,:,iTimeStep) = E;
    At_T(:,:,iTimeStep) = T;
    At_omega(:,:,iTimeStep) = omega;
    At_Kappa_vs(:,:,iTimeStep) = Kappa_vs;
    At_CTS(iTimeStep,:) = CTS;
    At_thk_TEMP(iTimeStep,:) = thk_TEMP;
    At_m_basal(iTimeStep,:) = m_basal;
    
    % calculate the basal water layer thickness
    if iTimeStep==1
        At_thk_w(iTimeStep,:) = zeros(1,M) + dt_u*m_basal;
    else
        At_thk_w(iTimeStep,:) = At_thk_w(iTimeStep-1,:) + dt_u*m_basal;
    end
    
    % logic value for basal thermal state
    logic1 = (At_thk_w(iTimeStep,:)>0 & At_isTEMP(iTimeStep,:)==0);
    At_isTEMP(iTimeStep,:) = At_isTEMP(iTimeStep,:) | logic1;
    
    if strcmpi(p.type_enth_solver,'MEGM')
        At_qw_TEMP_darcy(:,:,iTimeStep) = qw_TEMP_darcy;
        At_qw_TEMP(:,:,iTimeStep) = qw_TEMP;
    end   
end
toc

% figure
figure
subplot(1,3,1)
hold on
plot(E(:,1)/1000, H*zeta, 'k-', 'linewidth', 2)

xlabel('E (\times 10^3 J kg^{-1})')
ylabel('\zeta')
ylim([-5,205])
grid on
box on

subplot(1,3,2)
hold on
plot(T(:,1)-273.15, H*zeta, 'k-', 'linewidth', 2)

xlabel('T (^\circC)')
ylim([-5,205])
grid on
box on

subplot(1,3,3)
hold on
plot(omega(:,1)*100, H*zeta, 'k-', 'linewidth', 2)
xlabel('\omega (%)')
xlim([0, 4])
ylim([-5,205])
grid on
box on

set(gcf, 'position', [239 153 1005 412])