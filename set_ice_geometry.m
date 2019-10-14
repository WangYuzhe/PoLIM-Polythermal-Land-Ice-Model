function set_ice_geometry()

global xi dx hS hB H W W_s M Ms N dzeta zeta type_valley beta2_s

% Haut d'Arolla Glacier
arolla = dlmread('arolla100.dat');
arolla = arolla';

xi = arolla(1,:);
hB = arolla(2,:);
hS = arolla(3,:);
H = hS - hB + 0.1;

dx = xi(2) - xi(1);
M = length(xi);
Ms = M + 1;

N = 31;
dzeta = 1/(N-1);
zeta = 0:dzeta:1;
zeta = zeta';

beta2 = 1e4 + 1e4*sin(2*pi/5e3*xi);
beta2_s = main2staggerX(beta2);

% Slab Experiment
% L = 10*1000;
% M = 51;
% xi = linspace(0, L, M);
% slope = 1*pi/180; % mean slope [degree]
% hS = -xi.*tan(slope);
% hB = hS - 500;
% H = hS - hB;
% 
% dx = xi(2) - xi(1);
% Ms = M + 1;
% N = 31;
% dzeta = 1/(N-1);
% zeta = 0:dzeta:1;
% zeta = zeta';

% Half-width distribution
W = zeros(N,M);
switch type_valley
    case 1 % trapezoid
        Wsurf = linspace(1000,500,M);
%         W = ones(N,1)*Wsurf;
        alpha_trapezoid = 45*pi/180;
        for i = 1:M
            W(:,i) = Wsurf(i) - (H(i) - H(i)*zeta)/tan(alpha_trapezoid);
            
            % Make sure Wbase >= 1
            logic1 = (W(:,i) < 1e-3);
            W(:,i) = (1-logic1).*W(:,i) + 1*logic1;
        end

    case 2 % Svensson: y = ax^b
        Wsurf = 500*ones(1,M);
        for i=1:M
            W(:,i) = (zeta.^(1/q(i)))*Wsurf(i);
        end
        W(1,:) = W(2,:)/2;
    case 3 % rectangular
        Wsurf = 10000*ones(1,M);
        W = ones(N,1)*Wsurf;
end

W_s = [W(:,1), (W(:,2:end)+W(:,1:end-1))/2, W(:,end)];



end