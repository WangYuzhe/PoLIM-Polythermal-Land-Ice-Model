function [Hnp1] = get_evolution_continuity(Hn, u, smb, dt_u)
% The continuity equation is converted into a diffusion equation (Pattyn 2002,
% Pimentel 2011) for a rectangular cross section with variable width.
% Discretization uses a semi-implicit method and conservative boundary fluxes.

% Inputs:
% u: horizontal velocity field
% SMB: mass balance [m a-1]
% dt: time step [a]

% Outputs:
% Hnp1: ice thickness at the next time level

global M dzeta dx hS hB W dhSdx is_active i_term

if isempty(is_active)
    active = true(1, M);
    i_term_local = M;
else
    active = logical(is_active(:)');
    if numel(active) ~= M
        error('get_evolution_continuity:InvalidActiveMask', ...
            'is_active must contain one value per thickness grid point.')
    end
    i_term_local = find(active, 1, 'last');
    if isempty(i_term_local)
        Hnp1 = Hn;
        return
    end
    active(1:i_term_local) = true;
end
i_term = i_term_local;

Wsurf = W(end,:);
Wsurf_mid = (Wsurf(1:end-1) + Wsurf(2:end))/2;

uav = dzeta*(sum(u(1:end-1,:))+sum(u(2:end,:)))/2; % [1 by M]
uav_s = main2staggerX(uav);

diffu = zeros(1,M);
for i = 2:i_term-1
    diffu(i) = -uav_s(i+1).*Hn(i)./(dhSdx(i)+1e-10);
end
if i_term > 1
    diffu(1) = -uav_s(2).*Hn(1)./((hS(2)-hS(1))/dx+1e-10);
    diffu(i_term) = -uav_s(i_term+1).*Hn(i_term)./...
        (dhSdx(i_term)+1e-10);
end

diffu_mid = (diffu(1:end-1)+diffu(2:end))/2; % 1 by M-1

Hn_mid = (Hn(1:end-1) + Hn(2:end))/2;


LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);

% D = -uav * H / [(s(i+1) - s(i-1))/(2*dx)]
% semi-implicit, central difference
for i = 2:(i_term-1)
    LT1(i) = -diffu_mid(i-1)/(dx^2) +...
        diffu_mid(i)*(Wsurf(i+1)-Wsurf(i-1))/(4*dx^2*Wsurf(i));
    LT2(i) = 1/dt_u + (diffu_mid(i) + diffu_mid(i-1))/(dx^2);
    LT3(i) = -diffu_mid(i)/(dx^2) -...
        diffu_mid(i)*(Wsurf(i+1)-Wsurf(i-1))/(4*dx^2*Wsurf(i));
    RT(i) = Hn(i)/dt_u + diffu_mid(i)*(hB(i+1)-hB(i))/(dx^2) +...
        diffu_mid(i)*(hB(i+1)-hB(i-1))*...
        (Wsurf(i+1)-Wsurf(i-1))/(4*dx^2*Wsurf(i)) -...
        diffu_mid(i-1)*(hB(i)-hB(i-1))/(dx^2) + smb(i);
end

% Boundary nodes represent half control volumes; positive flux is downstream.
Q_head_in = uav_s(1)*Hn(1)*Wsurf(1);
Q_head_out = uav_s(2)*Hn_mid(1)*Wsurf_mid(1);
LT2(1) = 1;
RT(1) = Hn(1) + (Q_head_in-Q_head_out)*dt_u/...
    (Wsurf(1)*(dx/2)) + smb(1)*dt_u;

if i_term == 1
    RT(1) = Hn(1) + smb(1)*dt_u;
else
    Q_term_in = uav_s(i_term)*Hn_mid(i_term-1)*Wsurf_mid(i_term-1);
    Q_term_out = uav_s(i_term+1)*Hn(i_term)*Wsurf(i_term);
    LT2(i_term) = 1;
    RT(i_term) = Hn(i_term) + (Q_term_in-Q_term_out)*dt_u/...
        (Wsurf(i_term)*(dx/2)) + smb(i_term)*dt_u;
end

if i_term < M
    idx_inactive = (i_term+1):M;
    LT1(idx_inactive) = 0;
    LT2(idx_inactive) = 1;
    LT3(idx_inactive) = 0;
    RT(idx_inactive) = Hn(idx_inactive);
end

LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
Hnp1 = LT\RT; % [M by 1]
Hnp1 = Hnp1'; % [1 by M]
Hnp1(~active) = Hn(~active);
end
