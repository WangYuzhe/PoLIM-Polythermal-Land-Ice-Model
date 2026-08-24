function [hSnp1] = get_evolution_kinematic(hSn, u, w, SMB, dt)
% Relaxation of the free surface based on kinematic boundary.
% Date: 2016-8-25

global M dx

% dt: [a]

u_surf = u(end,:);
w_surf = w(end,:);

LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);

for i = 2:M-1
    LT1(i) = -u_surf(i)/(4*dx);
    LT2(i) = 1/dt;
    LT3(i) = u_surf(i)/(4*dx);
    RT(i) = u_surf(i)/(4*dx)*hSn(i-1) + hSn(i)/dt -...
        u_surf(i)/(4*dx)*hSn(i+1) + w_surf(i) + SMB(i);
end

LT1(1) = 0;
LT2(1) = 1/dt - u_surf(1)/(2*dx);
LT3(1) = u_surf(1)/(2*dx);
RT(1) = (1/dt + u_surf(1)/(2*dx))*hSn(1) -...
    u_surf(1)/(2*dx)*hSn(2) + w_surf(1) + SMB(1);

LT1(M) = -u_surf(end)/(2*dx);
LT2(M) = 1/dt + u_surf(end)/(2*dx);
LT3(M) = 0;
RT(M) = u_surf(end)/(2*dx)*hSn(M-1) +...
    (1/dt - u_surf(end)/(2*dx))*hSn(M) +...
    w_surf(end) + SMB(end);

LT = spdiags([[LT1(2:end);0],LT2,[0;LT3(1:end-1)]],[-1,0,1],M,M);
hSnp1 = LT\RT;

hSnp1 = hSnp1';
end