function a = compute_dart_drag_acceleration(r, v, et, params)
%COMPUTE_DART_DRAG_ACCELERATION Evaluate the DART ejecta-drag model.
% Naidu et al. (2024), Equation 4:
%   theta_ddot_drag = -A * exp(-(t - t_imp) / tau)   [rad/s^2]
%   a_drag = r * theta_ddot_drag * t_hat              [m/s^2]
% The negative tangential acceleration reduces the semimajor axis and period.
% The caller supplies the fitted impact epoch, amplitude, and decay time.

a = zeros(3,1);

if et < params.t_imp
    return;
end

dt_drag = et - params.t_imp;
if dt_drag > params.drag_cutoff
    return;
end

rn_cur = norm(r);
r_hat  = r / rn_cur;

% Normalize the transverse velocity.
v_tan   = v - dot(v, r_hat) * r_hat;
v_tan_n = norm(v_tan);
if v_tan_n <= 1e-12
    return;
end
t_hat = v_tan / v_tan_n;

% Convert the angular deceleration to tangential acceleration.
theta_ddot = -params.drag_A * exp(-dt_drag / params.drag_tau);
a = (rn_cur * theta_ddot) * t_hat;

end
