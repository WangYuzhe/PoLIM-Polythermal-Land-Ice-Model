function set_ice_geometry()

global xi dx hS hB H W W_s M Ms N dzeta zeta dzetadx
global type_valley type_geometry

switch type_geometry
    case 1 % Hault d'Arolla
        geoGlacier = dlmread('arolla100.dat');
        geoGlacier = geoGlacier';
        
        % origin
        xi = geoGlacier(1,:);
        hB = geoGlacier(2,:);
        hS = geoGlacier(3,:);
        H = hS - hB + 0.1;
        
        dx = xi(2) - xi(1);
        M = length(xi);
        Ms = M + 1;
        
        N = 31;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        
        % Half-width distribution
        W = zeros(N,M);
        switch type_valley
            case 1 % trapezoid
                alpha_trapezoid = 45*pi/180;
                Wsurf = 10*ones(1,M);
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
    case 2 % Hault d'Arolla (resampled to 250 points)
        geoGlacier = dlmread('arolla100.dat');
        geoGlacier = geoGlacier';
        
        % origin
        xi0 = geoGlacier(1,:);
        hB0 = geoGlacier(2,:);
        hS0 = geoGlacier(3,:);
        
        dx = 10;
        xi = 0:dx:5000;
        hB1 = interp1(xi0, hB0, xi, 'spline');
        hS1 = interp1(xi0, hS0, xi, 'spline');
        
        hB2 = smooth(hB1, 0.1, 'loess');
        hS2 = smooth(hS1, 0.1, 'loess');
        
        % plot(xi0, hB0, 'k', xi0, hS0, 'k')
        % hold on
        % % plot(xi, hB1, 'r', xi, hS1, 'r')
        % plot(xi, hB2', 'b', xi, hS2', 'b')
        
        % hS = hS1; hB = hB1;
        hS = hS2'; hB = hB2';
        
        logic1 = (hS-hB) > 0;
        hB = logic1.*hB + (1-logic1).*(hS - 0.1);
        H = hS - hB;
        
        M = length(xi);
        Ms = M + 1;
        
        N = 31;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        
        % Half-width distribution
        W = zeros(N,M);
        switch type_valley
            case 1 % trapezoid
                alpha_trapezoid = 45*pi/180;
                Wsurf = 10*ones(1,M);
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
    case 3 % Kleiner Exp. A (parallel-sided slab)
        H = 1000;         % ice thickness [m]
        N = 201;
        M = 1;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        dx = 100; % any constant not equal to zero
        dzetadx = zeros(N, M);
    case 4 % Kleiner Exp. B (polythermal parallel-sided slab)
        H = 200; % ice thickness [m]
        N = 401;  % N = 21, dz = 10 m; N = 401, dz = 0.5 m
        M = 1;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        dx = 100; % any constant not equal to zero
        dzetadx = zeros(N, M);
    case 5 % Hewitt-Schoof Ice Cap Experiment
        M = 51;
        Ms = M+1;
        L = 100*1000;
        xi = linspace(0, L, M);
        hB = zeros(1,M);
        hS = 1500.*(1 - (xi/L).^2);
        H = hS - hB + 0.1;
        dx = xi(2) - xi(1);
        
        N = 201; % dz = 200/(N-1), N = 200/dz + 1
        dzeta = 1/(N-1);
        zeta = linspace(0,1,N);
        zeta = zeta';
    case 6 % Storglaciaeren example
        geoGlacier = load('geoStorglaciaeren.mat');
        geoGlacier = geoGlacier.sg;
        geoGlacier = geoGlacier';
        
        % origin
        xi = geoGlacier(1, 5:104);
        hB = geoGlacier(2, 5:104);
        hS = geoGlacier(3, 5:104);
        H = geoGlacier(4, 5:104) + 0.1;
        
        dx = xi(2) - xi(1);
        M = length(xi);
        Ms = M + 1;
        
        N = 31;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        
        % Half-width distribution
        W = zeros(N,M);
        switch type_valley
            case 1 % trapezoid
                alpha_trapezoid = 45*pi/180;
                Wsurf = 10*ones(1,M);
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
                Wsurf = 1000*ones(1,M);
                W = ones(N,1)*Wsurf;
        end
        
        W_s = [W(:,1), (W(:,2:end)+W(:,1:end-1))/2, W(:,end)];
end
end