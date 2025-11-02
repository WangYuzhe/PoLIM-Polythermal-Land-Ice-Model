% Date: 2019-5-17
% Date: 2019-11-10

clc
clearvars
clearvars -global

global M xi rhow g hB dx

% parameters
p = params();
SPD = p.SPD;
SPY = p.SPY;

% ice geometry
ice_L = 100e3;
ice_W = 20e3;
set_ice_geometry('geo_shmip_icesheet', p);

% index for the three bands
Hindex1 = find(xi==10000); Hindex2 = find(xi==15000); % highest band
Mindex1 = find(xi==50000); Mindex2 = find(xi==55000); % middle band
Lindex1 = find(xi==85000); Lindex2 = find(xi==90000); % lower band

% TIME SETTING
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

% WATER SOURCES
m_basal = 7.93e-11*ones(1,M); % [m s-1]
m_moulin = zeros(1,M);
deltaT = [-4, -2, 0, 2, 4];

% PARAMETERS
ub = 1e-6*ones(1,M);
visc_b = 1e12*ones(1,M);
iter_max_period = 5;

for i = 1%:length(deltaT)
    fprintf('D%d is running!\n', i)
    m_runoff = calc_seasonRunoff(arrayTime, numTimeStep, deltaT(i));
    totalRechg = sum(m_runoff, 2)*ice_W*ice_L;
    
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
%                     load result_A1.mat
%                     Hw0 = Hw_A1;
%                     phi0 = phi_A1;
Hw0 = 1e-2*ones(1,M); % [m]
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
            [phi, Neff, Hw, qflux] = solver_subHydro1(dt, m_basal,...
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
    % Discharge for the band
    At_Qw = [mean(At_qflux(:, Hindex1:Hindex2),2)*ice_W,...
        mean(At_qflux(:, Mindex1:Mindex2),2)*ice_W,...
        mean(At_qflux(:, Lindex1:Lindex2),2)*ice_W]; % upper, middle, lower
    
    % SAVE RESULTS
    resultFileName = strcat('result_D', num2str(i), ' ');
    eval([ 'save ' resultFileName ' At_N ' 'At_qflux ' 'At_Qw ' 'totalRechg'])
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
