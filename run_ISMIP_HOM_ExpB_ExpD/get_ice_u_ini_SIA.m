% Date: 2017-8-19

function u_ini = get_ice_u_ini_SIA()
global AGlen rho g n dhSdx H N Ms zeta

u_ini = zeros(N, Ms);
for i = 1:Ms
    u_ini(:,i) = -2*AGlen*(rho*g*dhSdx(i))^n*((H(i) - H(i)*zeta).^(n+1)/(n+1) - H(i)^(n+1)/(n+1));    
end

end