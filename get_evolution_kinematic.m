% Relaxation of the free surface based on kinematic boundary.
% History:
% 2016-8-25

function [hS_new] = get_evolution_kinematic(u, w, B)
global hS hB H Hmin M dx

hS_lst = hS;
index = find(H < 0.2);

dt = 1; % [a]

us = u(end,:);
ws = w(end,:);

LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);

for i = 2:M-1
    LT1(i) = -us(i)/(4*dx);
    LT2(i) = 1/dt;
    LT3(i) = us(i)/(4*dx);
    RT(i) = us(i)/(4*dx)*hS_lst(i-1) + hS_lst(i)/dt -...
        us(i)/(4*dx)*hS_lst(i+1) + ws(i) + B(i);
end

LT1(1) = 0;
LT2(1) = 1/dt - us(1)/(2*dx);
LT3(1) = us(1)/(2*dx);
RT(1) = (1/dt + us(1)/(2*dx))*hS_lst(1) - us(1)/(2*dx)*hS_lst(2) + ws(1) + B(1);

LT1(M) = -us(end)/(2*dx);
LT2(M) = 1/dt + us(end)/(2*dx);
LT3(M) = 0;
RT(M) = us(end)/(2*dx)*hS_lst(M-1) + (1/dt - us(end)/(2*dx))*hS_lst(M) + ws(end) + B(end);

LT = spdiags([[LT1(2:end);0],LT2,[0;LT3(1:end-1)]],[-1,0,1],M,M);
hS_lst = LT\RT;

hS_lst(index) = hB(index) + Hmin;
    
hS_new = hS_lst';

end