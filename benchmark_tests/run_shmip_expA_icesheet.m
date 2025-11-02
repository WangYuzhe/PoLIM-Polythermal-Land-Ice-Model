clc
clearvars
clearvars -global

global M xi hB dx H

p = params();
SPD = p.SPD;
SPY = p.SPY;

% set_ice_geometry('./inputs/geo_hewitt_valley', p);
set_ice_geometry('./inputs/geo_shmip_glacier', p);

% TIME SETTING
dt_hydro = 0.5*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt_hydro, endTime);

% WATER SOURCES
% m_basal_range = [7.93e-11, 1.59e-9, 5.79e-9, 2.5e-8, 4.5e-8, 5.79e-7]; % [m s-1]
m_basal_range = 3.8e-10; % 12 mm a-1
m_moulin = zeros(1,M);
m_runoff = zeros(numTimeStep,M);
deltaT = [-4, -2, 0, 2, 4];

% PARAMETERS
ev = 0; % englacial void fraction
ub = 1e-6*ones(1,M);
visc_b = 1e12*ones(1,M);
iter_max_period = 20;

for i = 1%:length(m_basal_range)
    fprintf('D%d is running!\n', i)
    m_basal = m_basal_range(i)*ones(1,M);
    
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
                    Hw0 = 1e-4*ones(1,M); % [m]
                    phi0 = 1e5*ones(1,M); %rhow*g*hB;
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
            [phi, Neff, Hw, qflux] = solver_subHydro(dt_hydro, m_basal,...
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
    % SAVE RESULTS
    resultFileName = strcat('N_A', num2str(i));
    eval([resultFileName ' = Neff;'])
    eval(['save ' resultFileName ' ' resultFileName])
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
