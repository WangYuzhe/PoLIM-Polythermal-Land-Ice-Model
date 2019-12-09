function [Hn] = get_evolution_continuity(u, B)
% The continuity equation is converted into a diffusion equation (Pattyn 2002,
% Pimentel 2011).
% Discretization of the continuity equation uses the semi-implicit metHnd

% Inputs:
% u: horizontal velocity field
% B: mass balance

% Outputs:
% Hn: ice thickness

global M N dzeta dx hS hB H W dhSdx

type_formulation = 1;
Wsurf = W(end,:);
Hn = H;

dt = 1;
% dt_step = 1;
% dt0 = 0;

switch type_formulation
    case 1  % D = -uav * H / dhSdx
        % depth averaged horizontal velocity using trapezoidal rule
        uav = dzeta*(sum(u(1:N-1,:))+sum(u(2:N,:)))/2; % [1 by M]
        uav_s = [uav(1), (uav(1:end-1)+uav(2:end))/2, uav(end)];
        
        D = zeros(1,M);
        for i = 2:M-1
            D(i) = -uav_s(i+1).*Hn(i)./(dhSdx(i)+1e-10);
        end
        D(1) = -uav_s(2).*Hn(1)./((hS(2)-hS(1))/dx);
        D(M) = -uav_s(M+1).*Hn(M)./((hS(M)-hS(M-1))/dx);
        
        Dmid = (D(1:end-1)+D(2:end))/2; % 1 by M-1
        
        LT1 = zeros(M,1);
        LT2 = zeros(M,1);
        LT3 = zeros(M,1);
        RT = zeros(M,1);
        
        for i = 2:(M-1)
            LT1(i) = -Dmid(i)/(dx^2);
            LT2(i) = 1/dt + (3*Dmid(i) - Dmid(i-1))/(dx^2);
            LT3(i) = (-2*Dmid(i) + Dmid(i-1))/(dx^2);
            
            RT(i) = Hn(i)/dt + Dmid(i)*(hB(i+1)-2*hB(i)+hB(i-1))/(dx^2)...
                + (Dmid(i) - Dmid(i-1))*(hB(i+1)-hB(i))/(dx^2) + B(i);
        end
        LT2(1) = 1; RT(1) = Hn(1)/dt + B(1);
        LT2(M) = 1; RT(M) = Hn(M)/dt + B(M);
        
        LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
        
        Hn = LT\RT; % [M by 1]
        
    case 2 % D = -uav * H / [(H(i+1) - H(i-1))/(2*dx)]
        % depth averaged horizontal velocity using trapezoidal rule
        uav = dzeta*(sum(u(1:N-1,:))+sum(u(2:N,:)))/2;
        uav_s = [uav(1), (uav(1:end-1)+uav(2:end))/2, uav(end)];
                
        D = zeros(1,M); % 1 by M
        for i = 2:M-1
            D(i) = -uav_s(i).*Hn(i)./((Hn(i+1)-Hn(i-1))/(2*dx));
        end
        D(1) = -uav_s(1).*Hn(1)./((Hn(2)-Hn(1))/dx);
        D(M) = -uav_s(M).*Hn(M)./((Hn(M)-Hn(M-1))/dx);
        
        Dmid = (D(1:end-1)+D(2:end))/2; % 1 by M-1
        
        LT1 = zeros(M,1);
        LT2 = zeros(M,1);
        LT3 = zeros(M,1);
        RT = zeros(M,1);
        
        for i = 2:(M-1)
            LT1(i) = -Dmid(i)/(dx^2);
            LT2(i) = 1/dt + (3*Dmid(i) - Dmid(i-1))/(dx^2);
            LT3(i) = (-2*Dmid(i) + Dmid(i-1))/(dx^2);
            
            RT(i) = Hn(i)/dt + B(i);
        end
        LT2(1) = 1; RT(1) = Hn(1)/dt + B(1);
        LT2(M) = 1; RT(M) = Hn(M)/dt + B(M);
        
        LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
        
        Hn = LT\RT; % [M by 1]
        
    case 3  % D = -uav * H / [(s(i+1) - s(i-1))/(2*dx)]; width variation (forward difference)
        % depth averaged horizontal velocity using trapezoidal rule
        uav = dzeta*(sum(u(1:N-1,:))+sum(u(2:N,:)))/2; % [1 by M]
        uav_s = [uav(1), (uav(1:end-1)+uav(2:end))/2, uav(end)];
        
        D = zeros(1,M);
        for i = 2:M-1
            D(i) = -uav_s(i+1).*Hn(i)./(dhSdx(i)+1e-10);
        end
        D(1) = -uav_s(2).*Hn(1)./((hS(2)-hS(1))/dx);
        D(M) = -uav_s(M+1).*Hn(M)./((hS(M)-hS(M-1))/dx);
        
        Dmid = (D(1:end-1)+D(2:end))/2; % 1 by M-1
        
        LT1 = zeros(M,1);
        LT2 = zeros(M,1);
        LT3 = zeros(M,1);
        RT = zeros(M,1);
        
        for i = 2:(M-1)
            LT1(i) = -Dmid(i)/(dx^2);
            LT2(i) = 1/dt + (3*Dmid(i) - Dmid(i-1))/(dx^2) + Dmid(i)*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2);
            LT3(i) = (-2*Dmid(i) + Dmid(i-1))/(dx^2) - Dmid(i)*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2);
            
            RT(i) = Hn(i)/dt + Dmid(i)*(hB(i+1)-2*hB(i)+hB(i-1))/(dx^2)...
                + (Dmid(i) - Dmid(i-1))*(hB(i+1)-hB(i))/(dx^2) ...
                + Dmid(i)*(Wsurf(i+1)-Wsurf(i))*(hB(i+1)-hB(i))/Wsurf(i)/(dx^2) + B(i);
        end
        LT2(1) = 1; RT(1) = Hn(1)/dt + B(1);
        LT2(M) = 1; RT(M) = Hn(M)/dt + B(M);
        
        LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
        
        Hn = LT\RT; % [M by 1]
        
    case 4  % D = -uav * H / [(s(i+1) - s(i-1))/(2*dx)]; width variation (center difference)
        % depth averaged horizontal velocity using trapezoidal rule
        uav = dzeta*(sum(u(1:N-1,:))+sum(u(2:N,:)))/2; % [1 by M]
        uav_s = [uav(1), (uav(1:end-1)+uav(2:end))/2, uav(end)];
        
        D = zeros(1,M);
        for i = 2:M-1
            D(i) = -uav_s(i+1).*Hn(i)./(dhSdx(i)+1e-10);
        end
        D(1) = -uav_s(2).*Hn(1)./((hS(2)-hS(1))/dx);
        D(M) = -uav_s(M+1).*Hn(M)./((hS(M)-hS(M-1))/dx);
        
        Dmid = (D(1:end-1)+D(2:end))/2; % 1 by M-1
        
        LT1 = zeros(M,1);
        LT2 = zeros(M,1);
        LT3 = zeros(M,1);
        RT = zeros(M,1);
        
        alpha_comb = zeros(M,1);
        for i = 2:(M-1)
            alpha_comb(i) = Dmid(i)*(Wsurf(i+1)-Wsurf(i-1))/Wsurf(i)/(2*dx^2);
            LT1(i) = -Dmid(i)/(dx^2);
            LT2(i) = 1/dt + (3*Dmid(i) - Dmid(i-1))/(dx^2) + alpha_comb(i);
            LT3(i) = (-2*Dmid(i) + Dmid(i-1))/(dx^2) - alpha_comb(i);
            
            RT(i) = Hn(i)/dt + Dmid(i)*(hB(i+1)-2*hB(i)+hB(i-1))/(dx^2)...
                + (Dmid(i) - Dmid(i-1))*(hB(i+1)-hB(i))/(dx^2) ...
                + alpha_comb(i)*(hB(i+1)-hB(i)) + B(i);
        end
        LT2(1) = 1; RT(1) = Hn(1)/dt + B(1);
        LT2(M) = 1; RT(M) = Hn(M)/dt + B(M);
        
        LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
        
        Hn = LT\RT; % [M by 1]     
end

Hn = Hn'; % [1 by M]
end
