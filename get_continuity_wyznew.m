function [Hnp1] = get_continuity_wyznew(Hn, u, smb, dt_u)
% The continuity equation is converted into a diffusion equation (Pattyn 2002,
% Pimentel 2011).
% Discretization of the continuity equation uses the semi-implicit metHnd

% D = -uav * H / [(s(i+1) - s(i-1))/(2*dx)]
% semi-implicit, forward difference
% forward difference for width variation

% Inputs:
% u: horizontal velocity field
% SMB: mass balance [m a-1]
% dt: time step [a]

% Outputs:
% Hn: ice thickness

global M dzeta dx hS hB W dhSdx
eps_slp = -1e-6;

Wsurf = W(end,:);

% vertically averaged horizontal velocity [1 by M]
uav = dzeta*(sum(u(1:end-1,:))+sum(u(2:end,:)))/2;
uav_s = main2staggerX(uav);

diffu = zeros(1,M);
for i = 2:M-1
    diffu(i) = -uav_s(i+1).*Hn(i)./(dhSdx(i)+eps_slp);
end
diffu(1) = -uav_s(2).*Hn(1)./((hS(2)-hS(1))/dx+eps_slp);
diffu(end) = -uav_s(end).*Hn(end)./((hS(end)-hS(end-1))/dx+eps_slp);

diffu_mid = (diffu(1:end-1)+diffu(2:end))/2; % 1 by M-1

Hn_mid = (Hn(1:end-1) + Hn(2:end))/2;

LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);

for i = 2:M-1
    LT1(i) = -diffu_mid(i-1)/(dx^2);
    LT2(i) = 1/dt_u + (diffu_mid(i) + diffu_mid(i-1))/(dx^2) +...
        diffu_mid(i)*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2);
    LT3(i) = -diffu_mid(i)/(dx^2) -...
        diffu_mid(i)*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2);

    RT(i) = Hn(i)/dt_u + diffu_mid(i)*(hB(i+1)-hB(i))/(dx^2) -...
        diffu_mid(i-1)*(hB(i)-hB(i-1))/(dx^2) +...
        diffu_mid(i)*(hB(i+1)-hB(i))*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2) +...
        smb(i);
end

% boundary condition at glacier head
LT2(1) = 1; RT(1) = Hn(1) - uav_s(2)*Hn_mid(1)*dt_u/dx + smb(1)*dt_u;

% boundary condition at glacier terminus
LT2(end) = 1; RT(end) = Hn(end) + uav_s(end-1)*Hn_mid(end)*dt_u/dx + smb(end)*dt_u;

LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
Hnp1 = LT\RT; % [M by 1]
Hnp1 = Hnp1'; % [1 by M]
end
