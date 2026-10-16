function [q_B_to_J2000, omega_B, R_B_to_J2000, rtn] = ...
    initialize_dimorphos_attitude_3d(r_J2000, v_J2000, dot_phi_B)
%INITIALIZE_DIMORPHOS_ATTITUDE_3D Build the three-dimensional impact attitude.
% The synchronous pre-impact alignment is +X_B=-e_r, +Y_B=-e_theta,
% +Z_B=+h. The body frame subsequently evolves with the rigid body.

assert(isscalar(dot_phi_B) && isfinite(dot_phi_B), ...
    'cobird:attitude:initialize_dimorphos_attitude_3d:invalidLibrationRate', ...
    'dot_phi_B must be a finite scalar.');
[R_RTN_to_J2000, e_r, e_theta, h_hat, theta_dot] = ...
    cobird.attitude.build_rtn_frame(r_J2000, v_J2000);

R_B_to_J2000 = R_RTN_to_J2000 * diag([-1, -1, 1]);
q_B_to_J2000 = cobird.attitude.convert_dcm_to_quaternion(R_B_to_J2000);
omega_J2000 = (theta_dot + dot_phi_B) * h_hat;
omega_B = R_B_to_J2000' * omega_J2000;

rtn = struct('R_RTN_to_J2000', R_RTN_to_J2000, ...
    'e_r', e_r, 'e_theta', e_theta, 'h_hat', h_hat, ...
    'theta_dot', theta_dot);
end
