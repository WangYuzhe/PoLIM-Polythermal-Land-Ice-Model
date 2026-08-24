% Date: 2019-08-19
% Date: 2019-11-08
% Author: Wang Yuzhe
% using the expression for sheet evolution by Hoffman (2014)
% visc_b: viscosity of ice [Pa s], larger visc_ice_b, larger N; (old value: 1e12)

function [phi, Neff, Hw, qflux, rate_open, rate_close] = solver_subHydro(dt_hydro, m_basal,...
    m_moulin, m_runoff, Hw0, phi0, ub, visc_b, para)
global hB H dx M

%% PARAMETERS AND CONSTANTS
%
%
g = para.g;
rhoi = para.rhoi;
rhow = para.rhow;
ktransm = para.ktransm;
ev = para.ev;
eta_w = para.eta_w; % viscosity of water [Pa s]
hr = para.hr; % height of bedrock bumps [m]
lr = para.lr; % wavelength of bedrock bumps [m]

%% INITIALIZATIONS
%
%
grad_Hw0 = [(Hw0(2)-Hw0(1))/dx, (Hw0(3:end)-Hw0(1:end-2))/(2*dx),...
    (Hw0(end)-Hw0(end-1))/dx];
dHw3dx0 = 3*Hw0.^2.*grad_Hw0;

LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);
%% SOLUTION
%
%
for i = 2:M-1
    LT1(i) = -ktransm*Hw0(i)^3/(eta_w*dx^2) + ktransm*dHw3dx0(i)/(2*eta_w*dx);
    LT2(i) = 2*ktransm*Hw0(i)^3/(eta_w*dx^2) + Hw0(i)/visc_b(i) + ev/(rhow*g*dt_hydro);
    LT3(i) = -ktransm*Hw0(i)^3/(eta_w*dx^2) - ktransm*dHw3dx0(i)/(2*eta_w*dx);
    RT(i) = m_basal(i) + m_moulin(i) + m_runoff(i) -...
        ub(i)*(hr-Hw0(i))/lr +...
        Hw0(i)*(rhow*g*hB(i)+rhoi*g*H(i))/visc_b(i) +...
        ev*phi0(i)/(rhow*g*dt_hydro);
end

% upper boundary (zero flux condition)
LT1(1) = 0;
LT2(1) = 1;
LT3(1) = -1;
RT(1) = 0;

% lower boundary (Dirichlet condition: atmosphere pressure)
LT1(M) = 0;
LT2(M) = 1;
LT3(M) = 0;
RT(M) = 1e5; %rhoi*g*hB(end); % 1e5

% construct a sparse matrix
LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], M, M);

% solution
phi = LT\RT; % <M * 1>
phi = phi'; % <1 * M>

% phi(phi<rhoi*g*hB) = rhoi*g*hB(phi<rhoi*g*hB);

% effective pressure (0<=N<=Pi)
Pw = phi - rhoi*g*hB;
% Pw(Pw<0) = 0;
% phi(Pw<0) = rhoi*g*hB(Pw<0);

Pi = rhoi*g*H;
Neff = Pi - Pw;
Hw = Hw0 + (ub(i)*(hr-Hw0)/lr - Hw0.*Neff./visc_b)*dt_hydro;
% Neff(Neff<0) = 0;
% Neff(Neff>Pi) = Pi(Neff>Pi);

grad_phi = [(phi(2)-phi(1))/dx, (phi(3:end) - phi(1:end-2))/(2*dx), (phi(end)-phi(end-1))/dx];
qflux = -ktransm*Hw.^3.*grad_phi/eta_w; % [m2 s-1]

rate_open = ub.*(hr-Hw)./lr;
rate_close = Hw.*Neff./visc_b;


end