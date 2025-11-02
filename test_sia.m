
clear, clc

AGlen = 1e-16;
rhoi = 910;
g = 9.81;
n = 3;
slp = 4;
dhSdx = -tan(slp*pi/180);
H = 1.0;
zeta = linspace(0,1,31);
depth = H*(1-zeta);

sia = @(x) -2*AGlen/(n+1)*(rhoi*g*dhSdx)^n *...
    (H^(n+1) - (H - H*zeta(x)).^(n+1));

% zeta(N): surface; zeta(1): basal

N = length(zeta);

u = zeros(N,1);
for i=1:N
    u(i) = sia(i);
end

figure
plot(u, depth)
set(gca, 'YDir', 'reverse')

xlabel('velocity (m/a)')
ylabel('Depth (m)')

