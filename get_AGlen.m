function [AGlen_s] = get_AGlen(T)
% Date: 2017-12-17
% Author: Wang Yuzhe
% Calculate the flow rate factor A.

% Input:
% T: temperature

% Output:
% AGlen_s: flow rate factor A on staggered grid.

global N M H zeta SPY R betaCC type_Arrhenius

AGlen = zeros(N,M);
switch type_Arrhenius
    case 1 % Greve, 2009
        for i = 1:M
            A0_1 = 1.916e3; Q1 = 139000;
            A0_2 = 3.985e-13; Q2 = 60000;
            
            T_correct = T(:,i) + (1-zeta).*H(i)*betaCC;
            index1 = (T_correct > 263.15);
            index2 = (T_correct <= 263.15);
            
            AGlen(:,i) = A0_1*exp(-Q1./(R*T_correct)).*index1 +...
                A0_2*exp(-Q2./(R*T_correct)).*index2; % [Pa-3 s-1]
            AGlen(:,i) = AGlen(:,i)*SPY; % [Pa-3 a-1]
        end
    case 2 % Cuffey, 2010
        for i = 1:M
            T_correct = T(:,i) + (1-zeta).*H(i)*betaCC;
            index1 = (T_correct > 263);
            index2 = (T_correct <= 263);
            AGlen(:,i) = 5e-25*exp(-115e3/R*(1./T_correct - 1/263)).*index1 +...
                3.5e-25*exp(-60e3/R*(1./T_correct - 1/263)).*index2;
            AGlen(:,i) = AGlen(:,i)*SPY; % [Pa-3 a-1]
        end
end

AGlen_s = main2staggerX(AGlen);
end