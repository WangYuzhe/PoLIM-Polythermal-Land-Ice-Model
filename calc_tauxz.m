function [tauxz, tauxz_s] = calc_tauxz(u_s, visc_s)
% tau_xz = eta / H * du/dzeta
% calculate the basal shear stress (Pimentel & Flowers, 2010)
% Date: 2017-9-11

% Returns:
% tauxz (tauxz_s): basal shear stress [Pa]

% global H_s dzeta dx Ms M dzetadx_s
% ub = u(1,:);
% ub_s = u(1,:);
% 
% dudx_s = zeros(1,Ms);
% dudx_s(2:M) = (ub(2:end)-ub(1:end-1))/dx;
% dudx_s(1) = (ub_s(2)-ub_s(1))/dx;
% dudx_s(end) = (ub_s(end)-ub_s(end-1))/dx;
% 
% tauxz_s = visc_s(1,:)./H_s.*(-3*u_s(1,:)+4*u_s(2,:)-u_s(3,:))/(2*dzeta);
% 
% tauxz1_s = visc_s(1,:)./H_s.*(-3*u_s(1,:)+4*u_s(2,:)-u_s(3,:))/(2*dzeta) -...
%     4*visc_s(1,:).*(dudx_s + dzetadx_s(1,:).*(-3*u_s(1,:)+4*u_s(2,:)-u_s(3,:))/(2*dzeta));
% 
% tauxz = staggerX2main(tauxz_s);
% tauxz1 = staggerX2main(tauxz1_s);

global H_s dzeta

tauxz_s = visc_s(1,:)./H_s.*(-3*u_s(1,:)+4*u_s(2,:)-u_s(3,:))/(2*dzeta);
tauxz = staggerX2main(tauxz_s);

end