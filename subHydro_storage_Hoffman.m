% Date: 2019-08-19
% Date: 2019-11-08
% Author: Wang Yuzhe
% using the expression for sheet evolution by Hoffman (2014)
% visc_ice_b: viscosity of ice [Pa s], larger visc_ice_b, larger N; (old value: 1e12)

function [Neff, phi, qw, Hw, dHw3dx] = subHydro_storage_Hoffman(dt_hydro,...
    basalMelt, mouInput, seasonRunoff, Hw_nm1, dHw3dx_nm1, phi_nm1, ev, ub, visc_ice_b)
global rho rhow g hB H dx M

Ks = 4e-4; % transmissivity coefficient [], larger k0, larger N; (old value: 2e-4)
visc_water = 1e-3; % viscosity of water [Pa s]
hr = 0.1; % height of bedrock bumps [m]
lr = 2; % wavelength of bedrock bumps [m]

% INITIALIZATIONS
LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);

% RUN
for i = 2:M-1
    LT1(i) = -Ks*Hw_nm1(i)^3/(visc_water*dx^2) + Ks*dHw3dx_nm1(i)/(2*visc_water*dx);
    LT2(i) = 2*Ks*Hw_nm1(i)^3/(visc_water*dx^2) + Hw_nm1(i)/visc_ice_b(i) + ev/(rhow*g*dt_hydro);
    LT3(i) = -Ks*Hw_nm1(i)^3/(visc_water*dx^2) - Ks*dHw3dx_nm1(i)/(2*visc_water*dx);
    RT(i) = basalMelt(i) + mouInput(i) + seasonRunoff(i)...
        -ub(i)*(hr-Hw_nm1(i))/lr + Hw_nm1(i)/visc_ice_b(i)*(rhow*g*hB(i)+rho*g*H(i))...
        + ev*phi_nm1(i)/(rhow*g*dt_hydro);
end

% upper boundary (zero flux condition)
LT1(1) = 0; LT2(1) = 1; LT3(1) = -1; RT(1) = 0;

% lower boundary (Dirichlet condition: atmosphere pressure)
LT1(M) = 0; LT2(M) = 1; LT3(M) = 0; RT(M) = 1e5; %1e5

% construct a sparse matrix
LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], M, M);

% solution
phi = LT\RT; % <M * 1>
phi = phi'; % <1 * M>

% effective pressure
Neff = rhow*g*hB + rho*g*H - phi;
Hw = Hw_nm1 + (ub(i)*(hr-Hw_nm1)/lr - Hw_nm1.*Neff./visc_ice_b)*dt_hydro;
dHwdx = [(Hw(2)-Hw(1))/dx, (Hw(3:end) - Hw(1:end-2))/(2*dx), (Hw(end)-Hw(end-1))/dx];
dHw3dx = 3*Hw.^2.*dHwdx;

dphidx = [(phi(2)-phi(1))/dx, (phi(3:end) - phi(1:end-2))/(2*dx), (phi(end)-phi(end-1))/dx];
qw = -Ks*Hw.^3.*dphidx/visc_water; % [m2 s-1]

end