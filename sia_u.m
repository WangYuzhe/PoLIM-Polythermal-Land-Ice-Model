function [u_s] = sia_u(AGlen_s, para)
% solve the horizontal velocity field using shallow ice approximation (SIA)
% The solved velocity is used as an initial guess.
% Date: 2025/11/9

global dhSdx_s H_s zeta

rhoi = para.rhoi;
g = para.g;
n = para.n;

[N, Ms] = size(AGlen_s);
% sia = @(y,x) -2*AGlen_s(y,x)*(rhoi * g * dhSdx_s(x))^n/(n+1) *...
%     (H_s(x)^(n+1) - (H_s(x)-H_s(x)*zeta(y)).^(n+1));

u_s = zeros(N, Ms);
for i = 2:Ms-1
    local_slp = dhSdx_s(i);
    u_s(:,i) = -2*AGlen_s(:,i).*(rhoi * g * local_slp)^n/(n+1) .*...
        (H_s(i)^(n+1) - (H_s(i) - H_s(i)*zeta).^(n+1));
end

end