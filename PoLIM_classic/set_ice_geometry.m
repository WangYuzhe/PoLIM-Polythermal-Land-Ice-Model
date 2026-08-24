function set_ice_geometry(geofile, para)

global xi dx hS hB H W M Ms N dzeta zeta

type_valley = para.type_valley;

%% LOAD
%
%
load(geofile)
xi = geo.xi;
hB = geo.hB;
hS = geo.hS;
H = geo.H + para.Hmin;
% H = geo.H;
% H(end-1:end) = H(end-1:end) + para.Hmin;

%
% MODEL GRID
if size(xi)==1
    dx = 100;
else
    dx = xi(2) - xi(1);
end
M = length(xi);
Ms = M + 1;

N = para.layers;
dzeta = 1/(N-1);
zeta = 0:dzeta:1;
zeta = zeta';

%% WIDTH
%
Wsurf_const = 2e3;

% step 1: realistic plan-view half flowband width
% Wsurf: <1, M>
if isfield(geo, 'Wsurf')
    if size(geo.xi)==size(geo.Wsurf)
        Wsurf = geo.Wsurf;
        Wsurf(Wsurf<=0) = 0.5; % Handle the case where Wsurf=0 to avoid singularity (division by zero).
    else
        Wsurf = Wsurf_const/2*ones(1,M);
    end
else
    Wsurf = Wsurf_const/2*ones(1,M);
end

% step 2: half flowband width distribution
% for continuity, we use traditional variable names
W = zeros(N,M);
switch type_valley
    case 'sves' % in fact is 'parabola'. 
        W = sqrt(zeta)*Wsurf;
        W(1,:) = W(2,:)/2; % basal half-width
    case 'ellipse'
        W = sqrt(2*zeta - zeta.^2)*Wsurf;
        W(1,:) = W(2,:)/2; % basal half-width
    case 'trapz' % in fact is 'trapezoid'
        Wbasal = 5.0*ones(1,M); % basal half-width [m]
        W = ones(N,1)*Wbasal + zeta*(Wsurf - Wbasal);
    case 'rect' % in fact is 'rectangle'
        W = ones(N,1)*Wsurf;
    otherwise
        disp('type_valley should be sves, ellipse, trapz, rect!')
end

end