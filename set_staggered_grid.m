function set_staggered_grid()

global M Ms N dx zeta
global hS hB H hS_s hB_s H_s
global dhSdx_s dhSdx dhBdx_s dhBdx dzetadx_s dzetadx

hS_s = [hS(1), (hS(1:end-1)+hS(2:end))/2, hS(end)];
hB_s = [hB(1), (hB(1:end-1)+hB(2:end))/2, hB(end)];
H_s = [H(1), (H(1:end-1)+H(2:end))/2, H(end)];

dhSdx_s = zeros(1,Ms);
dhSdx_s(2:M) = (hS(2:end)-hS(1:end-1))/dx;
dhSdx = [dhSdx_s(1), (dhSdx_s(2:end-2)+dhSdx_s(3:end-1))/2, dhSdx_s(end)];

dhBdx_s = zeros(1,Ms);
dhBdx_s(2:M) = (hB(2:end)-hB(1:end-1))/dx;
dhBdx = [dhSdx_s(1), (dhBdx_s(2:end-2)+dhBdx_s(3:end-1))/2, dhBdx_s(end)];

dzetadx_s = zeros(N, Ms);  % dzetadx_s = f(H, hB, dzeta, dx)
for i = 1:Ms-1
    for j = 1:N
        
        if i == 1
            dzetadx_s(j,i) = -((-8*hB_s(i)+9*hB_s(i+1)-hB_s(i+2))/(3*dx) + zeta(j)*(-8*H_s(i)+9*H_s(i+1)-H_s(i+2))/(3*dx))/H_s(i);
        elseif i == 2
            dzetadx_s(j,i) = -((-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1))/(3*dx) + zeta(j)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1))/(3*dx))/H_s(i);
        elseif i == Ms-1
            dzetadx_s(j,i) = -((-1*hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1))/(3*dx) + zeta(j)*(-1*H_s(i-1)-3*H_s(i)+4*H_s(i+1))/(3*dx))/H_s(i);
        elseif i == Ms
            dzetadx_s(j,i) = -((8*hB_s(i)-9*hB_s(i-1)+hB_s(i-2))/(3*dx) + zeta(j)*(8*H_s(i)-9*H_s(i-1)+H_s(i-2))/(3*dx))/H_s(i);
        else
            dzetadx_s(j,i) = -((hB_s(i+1)-hB_s(i-1))/(2*dx) + zeta(j)*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
        end
        
    end
end
dzetadx = [dzetadx_s(:,1), (dzetadx_s(:,2:Ms-2)+dzetadx_s(:,3:Ms-1))/2, dzetadx_s(:,Ms)];

end