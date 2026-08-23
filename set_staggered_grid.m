function set_staggered_grid()

global M Ms N dx zeta
global hS hB H W hS_s hB_s H_s W_s
global dhSdx_s dhSdx dhBdx_s dhBdx dzetadx_s dzetadx
global is_active is_active_s i_term i_term_s

hS_s = main2staggerX(hS);
hB_s = main2staggerX(hB);
H_s = main2staggerX(H);
W_s = main2staggerX(W);

% Preserve the complete classic domain unless a driver supplies a fixed-grid
% activity mask. The mask is kept global to retain the classic architecture.
if isempty(is_active)
    is_active = true(1, M);
else
    is_active = logical(is_active(:)');
    if numel(is_active) ~= M
        error('set_staggered_grid:InvalidActiveMask', ...
            'is_active must contain one value per main-grid point.')
    end
end

i_term = find(is_active, 1, 'last');
if isempty(i_term)
    i_term = 0;
else
    % A locally inactive point must not split a physical flowline.
    is_active(1:i_term) = true;
end

is_active_s = false(1, Ms);
if i_term > 0
    i_term_s = min(i_term + 1, Ms);
    is_active_s(1:i_term_s) = true;
else
    i_term_s = 0;
end

dhSdx_s = zeros(1,Ms);
dhSdx_s(2:M) = (hS(2:end)-hS(1:end-1))/dx;
dhSdx_s(1) = (hS_s(2)-hS_s(1))/dx;
dhSdx_s(end) = (hS_s(end)-hS_s(end-1))/dx;

dhBdx_s = zeros(1,Ms);
dhBdx_s(2:M) = (hB(2:end)-hB(1:end-1))/dx;
dhBdx_s(1) = (hB_s(2)-hB_s(1))/dx;
dhBdx_s(end) = (hB_s(end)-hB_s(end-1))/dx;

dzetadx_s = zeros(N, Ms); % dzetadx_s = f(H, hB, dzeta, dx)
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

% Do not differentiate across the active--inactive jump. At a moving
% terminus, use the upstream geometry and zero derivatives on dormant cells.
if i_term_s > 0 && i_term < M
    if i_term > 1
        dhS_term = (hS(i_term) - hS(i_term - 1)) / dx;
        dhB_term = (hB(i_term) - hB(i_term - 1)) / dx;
        dH_term = (H(i_term) - H(i_term - 1)) / dx;
    else
        dhS_term = 0;
        dhB_term = 0;
        dH_term = 0;
    end

    dhSdx_s(i_term_s) = dhS_term;
    dhBdx_s(i_term_s) = dhB_term;
    dzetadx_s(:, i_term_s) = -(dhB_term + zeta*dH_term) ./ ...
        max(H_s(i_term_s), eps);

    dhSdx_s(i_term_s + 1:end) = 0;
    dhBdx_s(i_term_s + 1:end) = 0;
    dzetadx_s(:, i_term_s + 1:end) = 0;
elseif i_term_s == 0
    dhSdx_s(:) = 0;
    dhBdx_s(:) = 0;
    dzetadx_s(:) = 0;
end

dhSdx = staggerX2main(dhSdx_s);
dhBdx = staggerX2main(dhBdx_s);
dzetadx = staggerX2main(dzetadx_s);

dhSdx(~is_active) = 0;
dhBdx(~is_active) = 0;
dzetadx(:, ~is_active) = 0;

end
