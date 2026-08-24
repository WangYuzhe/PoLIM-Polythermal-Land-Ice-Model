function [hS, hB] = set_ice_geometry()

global experiment M Ms N dzeta zeta xi dx slope L
global dhSdx_s dhBdx_s dhBdx H H_s dzetadx_s dzetadx beta2_s

M = 41;
Ms = M + 1;
N = 41;
dzeta = 1/(N-1);
zeta = 0:dzeta:1;
zeta = zeta';

switch experiment
    case 1 % Experiment B
        L = 5*1000; % upper limit of x [km]: 160, 80, 40, 20, 10, 5
        xi = linspace(0, L, M);
        dx = xi(2)-xi(1);
        xi1 = [xi, xi(end)+dx, xi(end)+2*dx]; % 1, 2, ..., M, M+1, M+2 <---> 1, 2, ..., M, Ms, Ms+1
        
        slope = 0.5; % mean slope [degree]
        hS = -xi1.*tan(slope*pi/180); % <1 * M+2> = <1 * Ms+1>
        hB = hS - 1000 + 500*sin(2*pi/L*xi1); % <1 * M+2> = <1 * Ms+1>
        H = hS - hB; % <1 * M+2> = <1 * Ms+1>        
    case 2 % Experiment D
        L = 10*1000; % upper limit of x [km]: 160, 80, 40, 20, 10, 5
        xi = linspace(0, L, M);
        dx = xi(2)-xi(1);
        xi1 = [xi, xi(end)+dx, xi(end)+2*dx]; % 1, 2, ..., M, M+1, M+2 <---> 1, 2, ..., M, Ms, Ms+1
        
        slope = 0.1; % mean slope [degree]
        hS = -xi1.*tan(slope*pi/180);
        hB = hS - 1000;
        H = hS - hB;
        beta2 = 1000 + 1000*sin(2*pi/L*xi1); % [Pa a m-1]        
        beta2_s = (beta2(1:end-1)+beta2(2:end))/2;
end

dhSdx_s = (hS(2:end)-hS(1:end-1))/dx; % <1 * M+1> = <1 * Ms>
dhBdx_s = (hB(2:end)-hB(1:end-1))/dx; % <1 * M+1> = <1 * Ms>
H_s = (H(1:end-1)+H(2:end))/2; % <1 * M+1> = <1 * Ms>

dhBdx1 = (hB(3:end)-hB(1:end-2))/(2*dx);
dhBdx = [dhBdx1(end-1),dhBdx1];

dzetadx_s = zeros(N, M+1);
dzetadx = zeros(N, M+1);
for i = 1:Ms
    for j = 1:N
        dzetadx_s(j,i) = -((1-zeta(j))*dhBdx_s(i) + zeta(j)*dhSdx_s(i))/H_s(i);
        if i == 1
            dzetadx(j,i) = -((1-zeta(j))*(hB(M+1)-hB(M-1))/(2*dx) + zeta(j)*(hS(M+1)-hS(M-1))/(2*dx))/H(i);
        else
            dzetadx(j,i) = -((1-zeta(j))*(hB(i+1)-hB(i-1))/(2*dx) + zeta(j)*(hS(i+1)-hS(i-1))/(2*dx))/H(i);
        end
    end
end

end