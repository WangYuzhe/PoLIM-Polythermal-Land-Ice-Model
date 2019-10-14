clear, clc

SPY = 24*3600*365; % [s a-1]
AGlen = 2.4e-24;   % at 0 celcius [Pa-n s-1]
rho = 910;         % density of ice [kg m-3]
g = 9.81;          % accel of gravity [m s-2]
n = 3;             % exponent in Glen's flow law []

H = 500; % thickness [m]
slope = 1.1*pi/180; % surface slope

N = 41;
dzeta = 1/(N-1);
zeta = 0:dzeta:1;


u = -2*AGlen*(-rho*g*tan(slope)).^n * (-(H - H*zeta).^(n+1)/(n+1) + H^(n+1)/(n+1));
u = u*SPY;

plot(u, H*zeta, '-')
xlabel('Velocity (m a^{-1})')
ylabel('Elevation above bed (m)')