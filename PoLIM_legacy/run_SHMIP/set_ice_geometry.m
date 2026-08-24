function set_ice_geometry()

global xi dx hS hB H W W_s M Ms N dzeta zeta
global type_valley type_glacier_geometry

switch type_glacier_geometry
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
        
        Wsurf = 10000*ones(1,M);
        W = ones(N,1)*Wsurf;
        W_s = [W(:,1), (W(:,2:end)+W(:,1:end-1))/2, W(:,end)];
    case 2 % SHMIP ice sheet
        ice_L = 100e3;
        dx = 1000;
        xi = 0:dx:ice_L;
        M = length(xi);
        Ms = M + 1;
        
        hS = 6*((xi+5000).^(1/2)-sqrt(5000))+1;
        hB = zeros(1,M);
        hS = fliplr(hS);
        hB = fliplr(hB);
        H = hS - hB + 0.01;
        
        N = 31;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        Wsurf = 10000*ones(1,M);
        W = ones(N,1)*Wsurf;
        W_s = [W(:,1), (W(:,2:end)+W(:,1:end-1))/2, W(:,end)];        
    case 3 % SHMIP valley glacier
        ice_L = 6e3;
        dx = 60;
        xi = 0:dx:ice_L;
        M = length(xi);
        Ms = M + 1;
        
        gammab = 0.05;
        hS = 100*(xi+200).^(1/4) + xi/60 - (2*10^10)^(1/4) + 1;
        hB = (hS(end) - 6000*gammab) / (6000^2) * xi.^2 + gammab * xi;
        hS = fliplr(hS);
        hB = fliplr(hB);
        H = hS - hB + 1;
        
        N = 31;
        dzeta = 1/(N-1);
        zeta = 0:dzeta:1;
        zeta = zeta';
        Wsurf = 10000*ones(1,M);
        W = ones(N,1)*Wsurf;
        W_s = [W(:,1), (W(:,2:end)+W(:,1:end-1))/2, W(:,end)];        
end
end