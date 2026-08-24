function [hS, hB] = func_set_ice_geometry(L)

global experiment M Ms N dzeta zeta xi dx slope
global dhSdx_s dhBdx_s H_s beta2_s dzetadx_s dzetadx H dhBdx

M = 41;
Ms = M + 1;
N = 41;
dzeta = 1/(N-1);
zeta = 0:dzeta:1;
zeta = zeta';

switch experiment
    case 1 % Experiment B
        slope = 0.5; % mean slope [degree]
                
        xi = linspace(0, L, M);
        dx = xi(2)-xi(1);
        xi1 = [xi, xi(end)+dx, xi(end)+2*dx]; % 1, 2, ..., M, M+1, M+2 <---> 1, 2, ..., M, Ms, Ms+1
        hS = -xi1.*tan(slope*pi/180);
        hB = hS - 1000 + 500*sin(2*pi/L*xi1);
        H = hS - hB;
        
        dhSdx_s = (hS(2:end)-hS(1:end-1))/dx;
        dhBdx_s = (hB(2:end)-hB(1:end-1))/dx;
        H_s = (H(1:end-1)+H(2:end))/2;

        dhBdx1 = (hB(3:end)-hB(1:end-2))/(2*dx);
        dhBdx = [dhBdx1(end-1),dhBdx1];
        
        dzetadx_s = zeros(N, Ms);
        dzetadx = zeros(N, Ms);
        for i = 1:Ms
            for j = 1:N
                dzetadx_s(j,i) = -((1-zeta(j))*dhBdx_s(i) + zeta(j)*dhSdx_s(i))/H_s(i);
                if i == 1
                    dzetadx(j,i) = -((1-zeta(j))*(hB(Ms)-hB(Ms-2))/(2*dx) + zeta(j)*(hS(Ms)-hS(Ms-2))/(2*dx))/H(i);
                else
                    dzetadx(j,i) = -((1-zeta(j))*(hB(i+1)-hB(i-1))/(2*dx) + zeta(j)*(hS(i+1)-hS(i-1))/(2*dx))/H(i);
                end
            end
        end
    case 2 % Experiment D
        slope = 0.1; % mean slope [degree]
                
        xi = linspace(0, L, M);
        dx = xi(2)-xi(1);
        xi1 = [xi, xi(end)+dx, xi(end)+2*dx]; % 1, 2, ..., M, M+1, M+2 <---> 1, 2, ..., M, Ms, Ms+1
        hS = -xi1.*tan(slope*pi/180);
        hB = hS - 1000;
        H = hS - hB;
        beta2 = 1000 + 1000*sin(2*pi/L*xi1); % [Pa a m-1]
        
        dhSdx_s = (hS(2:end)-hS(1:end-1))/dx;
        dhBdx_s = (hB(2:end)-hB(1:end-1))/dx;
        H_s = (H(1:end-1)+H(2:end))/2;
        beta2_s = (beta2(1:end-1)+beta2(2:end))/2;
        
        dhBdx1 = (hB(3:end)-hB(1:end-2))/(2*dx);
        dhBdx = [dhBdx1(end-1),dhBdx1];
        
        dzetadx_s = zeros(N, Ms);
        dzetadx = zeros(N, Ms);
        for i = 1:Ms
            for j = 1:N
                dzetadx_s(j,i) = -((1-zeta(j))*dhBdx_s(i) + zeta(j)*dhSdx_s(i))/H_s(i);
                if i == 1
                    dzetadx(j,i) = -((1-zeta(j))*(hB(Ms)-hB(Ms-2))/(2*dx) + zeta(j)*(hS(Ms)-hS(Ms-2))/(2*dx))/H(i);
                else
                    dzetadx(j,i) = -((1-zeta(j))*(hB(i+1)-hB(i-1))/(2*dx) + zeta(j)*(hS(i+1)-hS(i-1))/(2*dx))/H(i);
                end
            end
        end
end

end