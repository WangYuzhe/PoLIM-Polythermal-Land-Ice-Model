function [E] = get_initial_enthalpy(Esbc)
% get the initial enthalpy field considering only thermal diffusion
% 2019-2-15

global kc Cp Qgeo M N dzeta H Kc

% initialization
LT1 = zeros(N,1);
LT2 = zeros(N,1);
LT3 = zeros(N,1);
RT = zeros(N,1);
E = zeros(N,M);

Kappa_s = Kc*ones(N-1,M);
for i = 1:M
        
    for j = 2:N-1
        RT(j) = 0;
        LT1(j) = -Kappa_s(j-1,i)/(H(i)^2*dzeta^2);
        LT2(j) = (Kappa_s(j-1,i)+Kappa_s(j,i))/(H(i)^2*dzeta^2);
        LT3(j) = -Kappa_s(j,i)/(H(i)^2*dzeta^2);
    end
    
    % surface BC
    LT1(N) = 0;
    LT2(N) = 1;
    LT3(N) = 0;
    RT(N) = Esbc(i);
    
    % basal BC
    LT1(1) = 0;
    LT2(1) = -1/dzeta;
    LT3(1) = 1/dzeta;
    RT(1) = -Qgeo*H(i)*Cp/kc;
    
    % left-hand system
    LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
    
    % solution
    E(:,i) = LT\RT; % <N * 1>
end

end
