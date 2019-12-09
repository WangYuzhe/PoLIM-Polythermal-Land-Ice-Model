function [Esbc] = set_thermalSBC()
% set the surface boundary condition for the thermal model
% 2019-2-15

global rho Cp Tref hS M

hS1 = 2500;
T1 = -2 + 273.15;

Tsbc = zeros(1,M);
for i = 1:M
    if hS(i) <= 2800
        Tsbc(i) = -0.006*(hS(i)-hS1) + T1;
    else
        Tsbc(i) = -2 + 273.15;
    end
end

% Esbc = rho*Cp*(Tsbc - Tref);

Esbc = Cp*(Tsbc - Tref);

end