% Date: 2017-9-11

function [Pdiff, Pdiff_s] = calc_Pdiff(u_s, visc_s)
global Ms dx dzetadx_s dzeta

Pdiff_s = zeros(1,Ms);

for i = 3:Ms-2
    Pdiff_s(i) = 2*visc_s(1,i)*(u_s(1,i+1)-u_s(1,i-1))/(2*dx) + ...
        2*visc_s(1,i)*dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta);
end

Pdiff_s(1) = 2*visc_s(1,1)*(u_s(1,Ms)-u_s(1,Ms-2))/(2*dx) + ...
    2*visc_s(1,1)*dzetadx_s(1,1)*(-3*u_s(1,1)+4*u_s(2,1)-u_s(3,1))/(2*dzeta);



Pdiff_s(Ms) = 2*visc_s(1,Ms)*(u_s(1,3)-u_s(1,1))/(2*dx) + ...
    2*visc_s(1,Ms)*dzetadx_s(1,Ms)*(-3*u_s(1,Ms)+4*u_s(2,Ms)-u_s(3,Ms))/(2*dzeta);

Pdiff_s = Pdiff_s/1000;

Pdiff = [Pdiff_s(:,1), (Pdiff_s(:, 2:Ms-2)+Pdiff_s(:, 3:Ms-1))/2, Pdiff_s(:,Ms)];

end