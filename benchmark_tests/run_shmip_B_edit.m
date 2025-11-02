clc
clearvars
clearvars -global

global M xi hB dx H

p = params();
SPD = p.SPD; % seconds in a day [s]
SPY = p.SPY; % seconds in a year [s]

set_ice_geometry('./inputs/geo_shmip_icesheet', p);

% TIME SETTING
dt_hydro = 0.5*SPD; % time step, 0.5/24*SPD [s]
startTime = 0;
endTime = 1*SPY; % [s]
[arrayTime, numTimeStep] = set_time_step(dt_hydro, startTime, endTime);

% WATER SOURCES
m_basal = 7.93e-11*ones(1,M); % basal melting
m_moulin = calc_expB_mouInput(4); % moulin input
m_runoff = zeros(numTimeStep,M); % distributed surface runoff input

% PARAMETERS
ev = 0; % englacial void fraction
ub = 1e-6*ones(1,M); % basal sliding velocity [m s-1]
visc_b = 1e12*ones(1,M); % basal ice viscosity [Pa s]
iter_max_period = 20;

% INITIAL CONDITIONS
At_Hw = zeros(numTimeStep, M);
At_phi = zeros(numTimeStep,M);
At_N = zeros(numTimeStep, M);
At_qflux = zeros(numTimeStep, M);

iter_period = 0;

while 1
    iter_period = iter_period + 1;
    for iTimeStep = 1:numTimeStep
        fprintf('Time Step: %d is running!\n', iTimeStep)
        
        if iter_period==1
            if iTimeStep==1
                Hw0 = 1e-3*ones(1,M); % initial sheet thickness [m]
                phi0 = 1e5*ones(1,M); % rhow*g*hB;
            else
                Hw0 = At_Hw(iTimeStep-1,:);
                phi0 = At_phi(iTimeStep-1,:);
            end
        else
            if iTimeStep==1
                Hw0 = At_Hw(end,:);
                phi0 = At_phi(end,:);
            else
                Hw0 = At_Hw(iTimeStep-1,:);
                phi0 = At_phi(iTimeStep-1,:);
            end
        end
        [phi, Neff, Hw, qflux] = solver_subHydro1(dt_hydro, m_basal,...
            m_moulin, m_runoff(iTimeStep,:), Hw0, phi0, ub, visc_b, p);
        
        At_Hw(iTimeStep, :) = Hw;
        At_phi(iTimeStep, :) = phi;
        At_N(iTimeStep, :) = Neff;
        At_qflux(iTimeStep, :) = qflux;
    end
    
    if iter_period>iter_max_period
        break
    end
end

figure
subplot(3,1,1)
hold on
F(1)=plot(xi/1e3, p.rhoi*p.g*H/1e6, 'k-');
F(2)=plot(xi/1e3, At_N(end,:)/1e6, 'r-');
F(3)=plot(xi/1e3, (phi-p.rhoi*p.g*hB)/1e6, 'b');
plot([xi(1)/1e3, xi(end)/1e3], [0,0], 'k--', 'Linewidth', 2)
xlabel('x (km)')
ylabel('Pressure (MPa)')
legend(F, 'Ice pressure', 'Effective pressure', 'Water pressure', 'Location', 'NorthEast')
box on

subplot(3,1,2)
plot(xi/1e3, At_Hw(end,:)*100, 'k-')
ylabel('Hw (cm)')
title('Sheet thickness')

subplot(3,1,3)
plot(xi/1e3, At_qflux(end,:), 'k-')
ylabel('q (m^2 s^{-1})')
title('Water flux')
set(gcf, 'position', [154.6000  111.1913  540.9391  683.6870])
