function [visc_s, visc] = get_ice_viscosity(u_s, u, AGlen_s)
global M Ms N dzeta dx de0 n dzetadx H

de2 = zeros(N, Ms); % [a-2]
visc = zeros(N, Ms); % [Pa a]

% At i = 1, j = N
de2(N,1) = ((3*u(N,1)-4*u(N-1,1)+u(N-2,1))/(2*dzeta))^2/(4*H(1)^2) + ...
    ((u_s(N,M)-u_s(N,M-1))/(dx) + dzetadx(N,1)*(3*u(N,1)-4*u(N-1,1)+u(N-2,1))/(2*dzeta))^2;
visc(N,1) = 1/2*AGlen_s(N,1)^(-1/n)*(de2(N,1) + de0)^((1-n)/(2*n));

% At i = 1, j = 1
de2(1,1) = ((-3*u(1,1)+4*u(2,1)-u(3,1))/(2*dzeta))^2/(4*H(1)^2) + ...
    ((u_s(N,M)-u_s(N,M-1))/(dx) + dzetadx(1,1)*(-3*u(1,1)+4*u(2,1)-u(3,1))/(2*dzeta))^2;
visc(1,1) = 1/2*AGlen_s(1,1)^(-1/n)*(de2(1,1) + de0)^((1-n)/(2*n));

% At i = Ms, j = N
de2(N,Ms) = ((3*u(N,Ms)-4*u(N-1,Ms)+u(N-2,Ms))/(2*dzeta))^2/(4*H(Ms)^2) + ...
    ((u_s(N,2)-u_s(N,1))/(dx) + dzetadx(N,Ms)*(3*u(N,Ms)-4*u(N-1,Ms)+u(N-2,Ms))/(2*dzeta))^2;
visc(N,Ms) = 1/2*AGlen_s(N,Ms)^(-1/n)*(de2(N,Ms) + de0)^((1-n)/(2*n));

% At i = Ms, j = 1
de2(1,Ms) = ((-3*u(1,Ms)+4*u(2,Ms)-u(3,Ms))/(2*dzeta))^2/(4*H(Ms)^2) + ...
    ((u_s(N,2)-u_s(N,1))/(dx) + dzetadx(1,Ms)*(-3*u(1,Ms)+4*u(2,Ms)-u(3,Ms))/(2*dzeta))^2;
visc(1,Ms) = 1/2*AGlen_s(1,Ms)^(-1/n)*(de2(1,Ms) + de0)^((1-n)/(2*n));

% At j = N & j = 1
for i = 2:Ms-1
    de2(N,i) = ((3*u(N,i)-4*u(N-1,i)+u(N-2,i))/(2*dzeta))^2/(4*H(i)^2) + ...
        ((u_s(N,i)-u_s(N,i-1))/(dx) + dzetadx(N,i)*(3*u(N,i)-4*u(N-1,i)+u(N-2,i))/(2*dzeta))^2;
    visc(N,i) = 1/2*AGlen_s(N,i)^(-1/n)*(de2(N,i) + de0)^((1-n)/(2*n));
    
    de2(1,i) = ((-3*u(1,i)+4*u(2,i)-u(3,i))/(2*dzeta))^2/(4*H(i)^2) + ...
        ((u_s(N,i)-u_s(N,i-1))/(dx) + dzetadx(1,i)*(-3*u(1,i)+4*u(2,i)-u(3,i))/(2*dzeta))^2;
    visc(1,i) = 1/2*AGlen_s(1,i)^(-1/n)*(de2(1,i) + de0)^((1-n)/(2*n));
end

% At i = 1 & i = Ms
for j = 2:N-1
    de2(j,1) = ((u(j+1,1)-u(j-1,1))/(2*dzeta))^2/(4*H(1)^2) + ...
        ((u_s(j,M)-u_s(j,M-1))/(dx) + dzetadx(j,1)*(u(j+1,1)-u(j-1,1))/(2*dzeta))^2;
    visc(j,1) = 1/2*AGlen_s(j,1)^(-1/n)*(de2(j,1) + de0)^((1-n)/(2*n));
    
    de2(j,Ms) = ((u(j+1,Ms)-u(j-1,Ms))/(2*dzeta))^2/(4*H(Ms)^2) + ...
        ((u_s(j,2)-u_s(j,1))/(dx) + dzetadx(j,Ms)*(u(j+1,Ms)-u(j-1,Ms))/(2*dzeta))^2;
    visc(j,Ms) = 1/2*AGlen_s(j,Ms)^(-1/n)*(de2(j,Ms) + de0)^((1-n)/(2*n));
end

% main
for i = 2:Ms-1
    for j = 2:N-1        
        de2(j,i) = ((u(j+1,i)-u(j-1,i))/(2*dzeta))^2/(4*H(i)^2) + ...
            ((u_s(j,i)-u_s(j,i-1))/(dx) + dzetadx(j,i)*(u(j+1,i)-u(j-1,i))/(2*dzeta))^2;
        visc(j,i) = 1/2*AGlen_s(j,i)^(-1/n)*(de2(j,i) + de0)^((1-n)/(2*n));
    end
end

visc1 = (visc(:,1:end-1)+visc(:,2:end))/2;
visc_s = [visc1, visc1(:,2)];
end