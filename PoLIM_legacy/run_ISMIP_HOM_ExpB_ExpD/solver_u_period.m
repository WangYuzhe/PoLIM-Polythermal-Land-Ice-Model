function [u_s] = solver_u_period(visc_s, visc)

global rho g AGlen n Ms N zeta dx dzeta experiment iter u_s_lst
global H_s dzetadx_s dhSdx_s dhBdx_s beta2_s dzetadx

L1 = zeros(N, Ms);
L2 = zeros(N, Ms);
L3 = zeros(N, Ms);
L4 = zeros(N, Ms);
L5 = zeros(N, Ms);

Ls = zeros(1, Ms);
Lb1 = zeros(1,Ms);
Lb2 = zeros(1,Ms);

LHS = zeros(N*Ms, N*Ms);
RHS = zeros(N*Ms, 1);

%------------------------------Build main----------------------------------
for i = 2:Ms-1
    L1(N,i) = 4*visc_s(N,i);
    L1(1,i) = 4*visc_s(1,i);
    
    L2(N,i) = 4*visc_s(N,i)*dzetadx_s(N,i)^2 + visc_s(N,i)/H_s(i)^2;
    L2(1,i) = 4*visc_s(1,i)*dzetadx_s(1,i)^2 + visc_s(1,i)/H_s(i)^2;
    
    L3(N,i) = 8*visc_s(N,i)*dzetadx_s(N,i);
    L3(1,i) = 8*visc_s(1,i)*dzetadx_s(1,i);
    
    L4(N,i) = 4*visc_s(N,i)*(dzetadx(N,i+1)-dzetadx(N,i))/(dx) +...
        4*visc_s(N,i)*dzetadx_s(N,i)*(3*dzetadx_s(N,i)-4*dzetadx_s(N-1,i)+dzetadx_s(N-2,i))/(2*dzeta) +...
        4*dzetadx_s(N,i)*(visc(N,i+1)-visc(N,i))/(dx) +...
        (4*dzetadx_s(N,i)^2+1/H_s(i)^2)*(3*visc_s(N,i)-4*visc_s(N-1,i)+visc_s(N-2,i))/(2*dzeta);
    L4(1,i) = 4*visc_s(1,i)*(dzetadx(1,i+1)-dzetadx(1,i))/(dx) +...
        4*visc_s(1,i)*dzetadx_s(1,i)*(-3*dzetadx_s(1,i)+4*dzetadx_s(2,i)-dzetadx_s(3,i))/(2*dzeta) +...
        4*dzetadx_s(1,i)*(visc(1,i+1)-visc(1,i))/(dx) +...
        (4*dzetadx_s(1,i)^2+1/H_s(i)^2)*(-3*visc_s(1,i)+4*visc_s(2,i)-visc_s(3,i))/(2*dzeta);
    
    L5(N,i) = 4*(visc(N,i+1)-visc(N,i-1))/(dx) +...
        4*dzetadx_s(N,i)*(3*visc_s(N,i)-4*visc_s(N-1,i)+visc_s(N-2,i))/(2*dzeta);
    L5(1,i) = 4*(visc(1,i+1)-visc(1,i))/(dx) +...
        4*dzetadx_s(1,i)*(-3*visc_s(1,i)+4*visc_s(2,i)-visc_s(3,i))/(2*dzeta);
    
    for j = 2:N-1
        L1(j,i) = 4*visc_s(j,i);
        L2(j,i) = 4*visc_s(j,i)*dzetadx_s(j,i)^2 + visc_s(j,i)/H_s(i)^2;
        L3(j,i) = 8*visc_s(j,i)*dzetadx_s(j,i);
        L4(j,i) = 4*visc_s(j,i)*(dzetadx(j,i+1)-dzetadx(j,i))/(dx) +...
            4*visc_s(j,i)*dzetadx_s(j,i)*(dzetadx_s(j+1,i)-dzetadx_s(j-1,i))/(2*dzeta) +...
            4*dzetadx_s(j,i)*(visc(j,i+1)-visc(j,i))/(dx) +...
            (4*dzetadx_s(j,i)^2+1/H_s(i)^2)*(visc_s(j+1,i)-visc_s(j-1,i))/(2*dzeta);
        L5(j,i) = 4*(visc(j,i+1)-visc(j,i))/(dx) +...
            4*dzetadx_s(j,i)*(visc_s(j+1,i)-visc_s(j-1,i))/(2*dzeta);
        
        k = (i-1)*N + j;
        RHS(k,1) = rho*g*dhSdx_s(i);
        
        LHS(k,k-N-1) = L3(j,i)/(4*dx*dzeta);
        LHS(k,k-N) = L1(j,i)/(dx^2) - L5(j,i)/(2*dx);
        LHS(k,k-N+1) = -L3(j,i)/(4*dx*dzeta);
        LHS(k,k-1) = L2(j,i)/(dzeta^2) - L4(j,i)/(2*dzeta);
        LHS(k,k) = -2*L1(j,i)/(dx^2) - 2*L2(j,i)/(dzeta^2);
        LHS(k,k+1) = L2(j,i)/(dzeta^2) + L4(j,i)/(2*dzeta);
        LHS(k,k+N-1) = -L3(j,i)/(4*dx*dzeta);
        LHS(k,k+N) = L1(j,i)/(dx^2) + L5(j,i)/(2*dx);
        LHS(k,k+N+1) = L3(j,i)/(4*dx*dzeta);
    end
end
%----------------------------Build surface BC (j=N)------------------------
for i = 1:Ms
    % Build RHS
    RHS(i*N,1) = rho*g*dhSdx_s(i);
    Ls(i) = 4*H_s(i)*dhSdx_s(i)/(1-4*H_s(i)*dzetadx_s(N,i)*dhSdx_s(i))*dzeta/dx;
end

% i == 2
LHS(2*N, N) = L1(N,2)/(dx^2) - Ls(2)*L2(N,2)/(dzeta^2) - Ls(2)*L4(N,2)/(2*dzeta) - L5(N,2)/(2*dx);
LHS(2*N, 2*N-1) = 2*L2(N,2)/(dzeta^2);
LHS(2*N, 2*N) = -2*L1(N,2)/(dx^2) - 2*L2(N,2)/(dzeta^2) - Ls(3)*L3(N,2)/(4*dx*dzeta);
LHS(2*N, 3*N) = L1(N,2)/(dx^2) + Ls(2)*L2(N,2)/(dzeta^2) + Ls(2)*L4(N,2)/(2*dzeta) + L5(N,2)/(2*dx);
LHS(2*N, 4*N) = Ls(3)*L3(N,2)/(4*dx*dzeta);
LHS(2*N, (Ms-2)*N) = Ls(Ms-1)*L3(N,2)/(4*dx*dzeta);
LHS(2*N, Ms*N) = -Ls(Ms-1)*L3(N,2)/(4*dx*dzeta);

% i == Ms-1
LHS((Ms-1)*N, N) = -Ls(2)*L3(N,Ms-1)/(4*dx*dzeta);
LHS((Ms-1)*N, 3*N) = Ls(2)*L3(N,Ms-1)/(4*dx*dzeta);
LHS((Ms-1)*N, (Ms-3)*N) = Ls(Ms-2)*L3(N,Ms-1)/(4*dx*dzeta);
LHS((Ms-1)*N, (Ms-2)*N) = L1(N,Ms-1)/(dx^2) - Ls(Ms-1)*L2(N,Ms-1)/(dzeta^2) -...
    Ls(Ms-1)*L4(N,Ms-1)/(2*dzeta) - L5(N,Ms-1)/(2*dx);
LHS((Ms-1)*N, (Ms-1)*N-1) = 2*L2(N,Ms-1)/(dzeta^2);
LHS((Ms-1)*N, (Ms-1)*N) = -2*L1(N,Ms-1)/(dx^2) - 2*L2(N,Ms-1)/(dzeta^2) -...
    Ls(Ms-2)*L3(N,Ms-1)/(4*dx*dzeta);
LHS((Ms-1)*N, Ms*N) = L1(N,Ms-1)/(dx^2) + Ls(Ms-1)*L2(N,Ms-1)/(dzeta^2) +...
    Ls(Ms-1)*L4(N,Ms-1)/(2*dzeta) + L5(N,Ms-1)/(2*dx);

% i = 3-->Ms-2
for i = 3:Ms-2
    LHS(i*N,(i-2)*N) = Ls(i-1)*L3(N,i)/(4*dx*dzeta);
    LHS(i*N,(i-1)*N) = L1(N,i)/(dx^2) - Ls(i)*L2(N,i)/(dzeta^2) - ...
        Ls(i)*L4(N,i)/(2*dzeta) - L5(N,i)/(2*dx);
    LHS(i*N,i*N-1) = 2*L2(N,i)/(dzeta^2);
    LHS(i*N,i*N) = -2*L1(N,i)/(dx^2) - 2*L2(N,i)/(dzeta^2) -...
        (Ls(i-1)+Ls(i+1))*L3(N,i)/(4*dx*dzeta);
    LHS(i*N,(i+1)*N) = L1(N,i)/(dx^2) + Ls(i)*L2(N,i)/(dzeta^2) +...
        Ls(i)*L4(N,i)/(2*dzeta) + L5(N,i)/(2*dx);
    LHS(i*N,(i+2)*N) = Ls(i+1)*L3(N,i)/(4*dx*dzeta);
end
%-----------------------------Build basal BC (j=1)-------------------------
switch experiment
    case 1 % No-slip (Exp. B)
        for i = 2:Ms-1
            LHS((i-1)*N+1,(i-1)*N+1) = 1;
            RHS((i-1)*N+1,1) = 0;
        end
        
    case 2 % Sliding (Exp. D)
        for i = 1:Ms
            RHS((i-1)*N+1,1) = rho*g*dhSdx_s(i);
            Lb1(i) = 4*dzeta*H_s(i)*dhBdx_s(i)/dx/(1-4*H_s(i)*dzetadx_s(1,i)*dhBdx_s(i));
            Lb2(i) = 2*dzeta*H_s(i)*beta2_s(i)/visc_s(1,i)/(1-4*H_s(i)*dzetadx_s(1,i)*dhBdx_s(i));
        end
        
        % i = 2
        LHS(N+1, 1) = L1(1,2)/(dx^2) + Lb1(2)*L2(1,2)/(dzeta^2) - ...
            Lb1(2)*L4(1,2)/(2*dzeta) - L5(1,2)/(2*dx);
        LHS(N+1, N+1) = -2*L1(1,2)/(dx^2) - (2+Lb2(2))*L2(1,2)/(dzeta^2) - ...
            Lb1(3)*L3(1,2)/(4*dx*dzeta) + Lb2(2)*L4(1,2)/(2*dzeta);
        LHS(N+1, N+2) = 2*L2(1,2)/(dzeta^2);
        LHS(N+1, 2*N+1) = L1(1,2)/(dx^2) - Lb1(2)*L2(1,2)/(dzeta^2) + ...
            Lb2(3)*L3(1,2)/(4*dx*dzeta) + Lb1(2)*L4(1,2)/(2*dzeta) + L5(1,2)/(2*dx);
        LHS(N+1, 3*N+1) = Lb1(3)*L3(1,2)/(4*dx*dzeta);
        LHS(N+1, (Ms-3)*N+1) = Lb1(Ms-1)*L3(1,2)/(4*dx*dzeta);
        LHS(N+1, (Ms-2)*N+1) = -Lb2(Ms-1)*L3(1,2)/(4*dx*dzeta);
        LHS(N+1, (Ms-1)*N+1) = -Lb1(Ms-1)*L3(1,2)/(4*dx*dzeta);
        
        % i = Ms-1
        LHS((Ms-2)*N+1, 1) = -Lb1(2)*L3(1,Ms-1)/(4*dx*dzeta);
        LHS((Ms-2)*N+1, N+1) = Lb2(2)*L3(1,Ms-1)/(4*dx*dzeta);
        LHS((Ms-2)*N+1, 2*N+1) = Lb1(2)*L3(1,Ms-1)/(4*dx*dzeta);
        LHS((Ms-2)*N+1, (Ms-4)*N+1) = Lb1(Ms-2)*L3(1,Ms-1)/(4*dx*dzeta);
        LHS((Ms-2)*N+1, (Ms-3)*N+1) = L1(1,Ms-1)/(dx^2) + Lb1(Ms-1)*L2(1,Ms-1)/(dzeta^2) - ...
            Lb2(Ms-2)*L3(1,Ms-1)/(4*dx*dzeta) - Lb1(Ms-1)*L4(1,Ms-1)/(2*dzeta) - L5(1,Ms-1)/(2*dx);
        LHS((Ms-2)*N+1, (Ms-2)*N+1) = -2*L1(1,Ms-1)/(dx^2) - (Lb2(Ms-1)+2)*L2(1,Ms-1)/(dzeta^2) -...
            Lb1(Ms-2)*L3(1,Ms-1)/(4*dx*dzeta) + Lb2(Ms-1)*L4(1,Ms-1)/(2*dzeta);
        LHS((Ms-2)*N+1, (Ms-2)*N+2) = 2*L2(1,Ms-1)/(dzeta^2);
        LHS((Ms-2)*N+1, (Ms-1)*N+1) = L1(1,Ms-1)/(dx^2) - Lb1(Ms-1)*L2(1,Ms-1)/(dzeta^2) + ...
            Lb1(Ms-1)*L4(1,Ms-1)/(2*dzeta) + L5(1,Ms-1)/(2*dx);
        
        % i = 3-->Ms-2
        for i = 3:Ms-2
            LHS((i-1)*N+1, (i-3)*N+1) = Lb1(i-1)*L3(1,i)/(4*dx*dzeta);
            LHS((i-1)*N+1, (i-2)*N+1) = L1(1,i)/(dx^2) + Lb1(i)*L2(1,i)/(dzeta^2) - ...
                Lb2(i-1)*L3(1,i)/(4*dx*dzeta) - Lb1(i)*L4(1,i)/(2*dzeta) - L5(1,i)/(2*dx);
            LHS((i-1)*N+1, (i-1)*N+1) = -2*L1(1,i)/(dx^2) - (2+Lb2(i))*L2(1,i)/(dzeta^2) -...
                (Lb1(i+1)+Lb1(i-1))*L3(1,i)/(4*dx*dzeta) + Lb2(i)*L4(1,i)/(2*dzeta);
            LHS((i-1)*N+1, (i-1)*N+2) = 2*L2(1,i)/(dzeta^2);
            LHS((i-1)*N+1, i*N+1) = L1(1,i)/(dx^2) - Lb1(i)*L2(1,i)/(dzeta^2) +...
                Lb2(i+1)*L3(1,i)/(4*dx*dzeta) + Lb1(i)*L4(1,i)/(2*dzeta) + L5(1,i)/(2*dx);
            LHS((i-1)*N+1, (i+1)*N+1) = Lb1(i+1)*L3(1,i)/(4*dx*dzeta);
        end
        
end
%------------------------Build lateral BC (i=1, i=Ms)---------------------
% periodic boundary condition
for j = 1:N
    if iter == 1 % initial guess
        % i = 1, upper lateral
        LHS(j,j) = 1;
%         RHS(j,1) = -2*AGlen*(rho*g*dhSdx_s(1))^n*((H_s(1) - H_s(1)*zeta(j)).^(n+1)/(n+1) - H_s(1)^(n+1)/(n+1));
        RHS(j,1) = 11;
        % i = Ms, lower lateral
        LHS((Ms-1)*N+j,(Ms-1)*N+j) = 1;
%         RHS((Ms-1)*N+j,1) = -2*AGlen*(rho*g*dhSdx_s(Ms))^n*((H_s(Ms) - H_s(Ms)*zeta(j)).^(n+1)/(n+1) - H_s(Ms)^(n+1)/(n+1));
        RHS((Ms-1)*N+j,1) = 11;
    else
        % i = 1 <---> i = Ms-1
        LHS(j,j) = 1;
        RHS(j,1) = u_s_lst(j,Ms-1);
        
        % i = Ms <---> i = 2
        LHS((Ms-1)*N+j,(Ms-1)*N+j) = 1;
        RHS((Ms-1)*N+j,1) = u_s_lst(j,2);
    end    
end
%---------------------------Solve linear system----------------------------
v = LHS\RHS;

u_s = zeros(N, Ms);
for k = 1:N*Ms    
    if fix(k/N) == k/N
        i = k/N;
    else
        i = fix(k/N) + 1;
    end    
    j = k-(i-1)*N;
    u_s(j,i) = v(k);
end

end