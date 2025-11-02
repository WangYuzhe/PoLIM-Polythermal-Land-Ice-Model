clear,clc

% Hmin = 0.01; % minimum thickness in case numerical singularity
% geo.xi = linspace(0, 750e3, 21);
% M = length(geo.xi);
% geo.hB = zeros(1, M);
% geo.hS = geo.hB + Hmin;
% geo.H = geo.hS - geo.hB;
% 
% save geo_eismint1 geo



% Parameters
n = 3;                % Glen's law exponent
R = 75e3;            % ice cap radius [m]
H0 = 100.0;            % central ice thickness [m]
M = 51;             % number of grid points

% Grid (0 at dome center, R at margin)
xi = linspace(0, R, M);

% Exponents
p = (n+1)/n;          % = 4/3 for n=3
q = n/(2*n+2);        % = 3/8 for n=3

% Halfar dome thickness profile
H = H0 * (1 - (xi/R).^p).^q;

% Bedrock elevation
hB = zeros(1,M);

% Surface elevation (bed = 0 here)
hS = hB + H;  

% Save it to mat
geo.xi = xi;
geo.hS = hS;
geo.hB = hB;
geo.H = H;

save geo_eismint geo
