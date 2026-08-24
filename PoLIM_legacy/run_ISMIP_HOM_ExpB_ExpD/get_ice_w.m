function [w_s, w] = get_ice_w(u_s, u)

global Ms N H dx dzeta dzetadx dhBdx

%----------------------------------U integrand-----------------------------
U = zeros(N,Ms);
for i = 2:Ms-1
    U(N,i) = (u_s(N,i)-u_s(N,i-1))/(dx) + dzetadx(N,i)*(3*u(N,i)-4*u(N-1,i)+u(N-2,i))/(2*dzeta);
    U(1,i) = (u_s(1,i)-u_s(1,i-1))/(dx) + dzetadx(1,i)*(-3*u(1,i)+4*u(2,i)-u(3,i))/(2*dzeta);
end

for j = 2:N-1
    U(j,1) = (u_s(j,Ms-1)-u_s(j,Ms-2))/(dx) + dzetadx(j,1)*(u(j+1,1)-u(j-1,1))/(2*dzeta);
    U(j,Ms) = (u_s(j,2)-u_s(j,1))/(dx) + dzetadx(j,Ms)*(u(j+1,Ms)-u(j-1,Ms))/(2*dzeta);
end

U(1,1) = (u_s(1,Ms-1)-u_s(1,Ms-2))/(dx) + dzetadx(1,1)*(-3*u(1,1)+4*u(2,1)-u(3,1))/(2*dzeta);
U(1,Ms) = (u_s(1,2)-u_s(1,1))/(dx) + dzetadx(1,Ms)*(-3*u(1,Ms)+4*u(2,Ms)-u(3,Ms))/(2*dzeta);
U(N,1) = (u_s(N,Ms-1)-u_s(N,Ms-2))/(dx) + dzetadx(N,1)*(3*u(N,1)-4*u(N-1,1)+u(N-2,1))/(2*dzeta);
U(N,Ms) = (u_s(N,2)-u_s(N,1))/(dx) + dzetadx(N,Ms)*(3*u(N,Ms)-4*u(N-1,Ms)+u(N-2,Ms))/(2*dzeta);

for i = 2:Ms-1
    for j = 2:N-1
        U(j,i) = (u_s(j,i)-u_s(j,i-1))/(dx) + dzetadx(j,i)*(u(j+1,i)-u(j-1,i))/(2*dzeta);
    end
end
%-----------------------------vertical velocity----------------------------
w_s = zeros(N,Ms);
for i = 1:Ms
    w_s(1,i) = u(1,i)*dhBdx(i) - H(i)*(U(1,i)/2 + 0*sum(U(2:j,i)))*dzeta;
    for j = 2:N
        w_s(j,i) = u(1,i)*dhBdx(i) - H(i)*(U(1,i)/2 + sum(U(2:j,i)))*dzeta;
    end
end

w = zeros(N,Ms);
for i = 1:Ms
    w(:,i) = [u(1,i)*dhBdx(i); (w_s(1:end-1,i)+w_s(2:end,i))/2]; % wb = ub*[partial(b)/partial(x)], ub=u(1,i)
end

end