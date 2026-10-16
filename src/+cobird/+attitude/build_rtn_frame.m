function [R_RTN_to_J2000, e_r, e_theta, h_hat, theta_dot] = build_rtn_frame(r_J2000, v_J2000)
%BUILD_RTN_FRAME Build a right-handed RTN frame from relative state vectors.
% Columns point along Didymos-to-Dimorphos (+R), prograde motion (+T),
% and orbital angular momentum (+N).

r_J2000 = r_J2000(:);
v_J2000 = v_J2000(:);
assert(numel(r_J2000) == 3 && numel(v_J2000) == 3, ...
    'cobird:attitude:build_rtn_frame:invalidStateDimensions', ...
    'Position and velocity must be three-element vectors.');

r_norm = norm(r_J2000);
h_vec = cross(r_J2000, v_J2000);
h_norm = norm(h_vec);
assert(isfinite(r_norm) && r_norm > 0, ...
    'cobird:attitude:build_rtn_frame:invalidPosition', ...
    'The relative position vector must have a finite, positive norm.');
assert(isfinite(h_norm) && h_norm > 0, ...
    'cobird:attitude:build_rtn_frame:invalidAngularMomentum', ...
    'The orbital angular-momentum vector must have a finite, positive norm.');

e_r = r_J2000 / r_norm;
h_hat = h_vec / h_norm;
e_theta = cross(h_hat, e_r);
e_theta = e_theta / norm(e_theta);
R_RTN_to_J2000 = [e_r, e_theta, h_hat];
theta_dot = h_norm / r_norm^2;
end
