function [dt_CFL, umax] = calc_thk_CFL_step(u_s, p, dx)
%CALC_THK_CFL_STEP Calculate a stable thickness-evolution time step.
%
% The continuity solver transports ice with depth-averaged staggered
% velocity, so its CFL condition must use the same velocity and only the
% faces belonging to the active glacier domain.

global N Ms zeta i_term

if ~isscalar(dx) || ~isfinite(dx) || dx <= 0
    error('calc_thk_CFL_step:InvalidGridSpacing', ...
        'dx must be a positive finite scalar.')
end

required = {'CFL', 'dt_H_min', 'dt_H_max'};
has_valid_params = all(isfield(p, required));
if has_valid_params
    vals = [p.CFL, p.dt_H_min, p.dt_H_max];
    has_valid_params = numel(vals) == 3 && all(isfinite(vals)) && ...
        all(vals > 0) && p.dt_H_max >= p.dt_H_min;
end
if ~has_valid_params
    error('calc_thk_CFL_step:InvalidParameters', ...
        ['CFL, dt_H_min, and dt_H_max must be positive scalars, and ' ...
         'dt_H_max must not be smaller than dt_H_min.'])
end

if ~isequal(size(u_s), [N, Ms]) || any(~isfinite(u_s(:)))
    error('calc_thk_CFL_step:InvalidVelocity', ...
        'u_s must be a finite N-by-Ms staggered velocity field.')
end

if isempty(i_term) || i_term == 0
    umax = 0;
    dt_CFL = p.dt_H_max;
    return
end

uav_s = trapz(zeta, u_s, 1);
idx_faces = 1:min(i_term + 1, Ms);
umax = max(abs(uav_s(idx_faces)));

if umax > 0
    dt_CFL = p.CFL*dx/umax;
else
    dt_CFL = p.dt_H_max;
end

% Raising the step to dt_H_min would violate the CFL condition.
if dt_CFL < p.dt_H_min
    error('calc_thk_CFL_step:StepBelowMinimum', ...
        ['Required CFL step %.6g yr is below dt_H_min %.6g yr ' ...
         '(maximum depth-averaged speed %.6g m/yr, dx %.6g m).'], ...
        dt_CFL, p.dt_H_min, umax, dx)
end

end
