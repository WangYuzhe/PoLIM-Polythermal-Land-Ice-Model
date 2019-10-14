function [w_vs, w] = get_ice_w(u_s, u)

global M N dx dzeta H dzetadx dhBdx

U = zeros(N,M);
    % s-1
w_vs = zeros(N,M);
w = zeros(N,M);
%------------------------------Build U integrand---------------------------
% You should know the essence of staggered grid in xi axis.
for i = 2:M-1
    for j = 2:N-1
        U(j,i) = (u_s(j,i+1)-u_s(j,i))/dx + dzetadx(j,i)*(u(j+1,i)-u(j-1,i))/(2*dzeta);
        U(N,i) = (u_s(N,i+1)-u_s(N,i))/dx + dzetadx(N,i)*(u(N,i)-u(N-1,i))/(dzeta);
        U(1,i) = (u_s(1,i+1)-u_s(1,i))/dx + dzetadx(1,i)*(u(2,i)-u(1,i))/(dzeta);
        U(j,1) = (u_s(j,2)-u_s(j,1))/(dx/2) + dzetadx(j,1)*(u(j+1,1)-u(j-1,1))/(2*dzeta);
        U(j,M) = (u_s(j,M+1)-u_s(j,M))/(dx/2) + dzetadx(j,M)*(u(j+1,M)-u(j-1,M))/(2*dzeta);
    end
end
U(1,1) = (u_s(1,2)-u_s(1,1))/(dx/2) + dzetadx(1,1)*(u(2,1)-u(1,1))/(dzeta);
U(1,M) = (u_s(1,M+1)-u_s(1,M))/(dx/2) + dzetadx(1,M)*(u(2,M)-u(1,M))/(dzeta);
U(N,1) = (u_s(N,2)-u_s(N,1))/(dx/2) + dzetadx(N,1)*(u(N,1)-u(N-1,1))/(dzeta);
U(N,M) = (u_s(N,M+1)-u_s(N,M))/(dx/2) + dzetadx(N,M)*(u(N,M)-u(N-1,M))/(dzeta);
%
%---------------------Calculate the vertical velocity----------------------
for i = 2:M-1
    for j = 2:N
        w_vs(j,i) = u(1,i)*dhBdx(i) - H(i)*(U(1,i)/2 + sum(U(2:j,i)))*dzeta;    % Covering grids, (j=2:N)~(i=2:M-1)
        w_vs(1,i) = u(1,i)*dhBdx(i) - H(i)*(U(1,i)/2 + 0*sum(U(1:2,i)))*dzeta;  % Base(j=1), i=2:M-1
        w_vs(j,1) = u(1,1)*dhBdx(1) - H(1)*(U(1,1)/2 + sum(U(2:j,1)))*dzeta;    % Head(i=1), j=2:N
        w_vs(j,M) = u(1,M)*dhBdx(M) - H(M)*(U(1,M)/2 + sum(U(2:j,M)))*dzeta;    % Terminus(i=M), j=2:N
    end
end
w_vs(1,1) = u(1,1)*dhBdx(1) - H(1)*(U(1,1)/2 + 0*sum(U(1:2,1)))*dzeta;    % Lower left(j=1,i=1)
w_vs(1,M) = u(1,M)*dhBdx(M) - H(M)*(U(1,M)/2 + 0*sum(U(1:2,M)))*dzeta;    % Lower right(j=1, i=M)
%
%--------------------------------------------------------------------------
for i = 2:M-1
    w(:,i) = [u(1,i)*dhBdx(i); (w_vs(1:end-1,i)+w_vs(2:end,i))/2];    % wb = ub*[partial(b)/partial(x)], ub=u(1,i)
end
w(:,1) = [u(1,1)*dhBdx(1); (w_vs(1:end-1,1)+w_vs(2:end,1))/2];
w(:,M) = [u(1,M)*dhBdx(M); (w_vs(1:end-1,M)+w_vs(2:end,M))/2];

end