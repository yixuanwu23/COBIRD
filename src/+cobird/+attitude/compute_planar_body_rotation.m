function rotation_body_to_inertial = compute_planar_body_rotation(theta, phi_B, e_r0, e_th0, h0_hat)
%COMPUTE_PLANAR_BODY_ROTATION Return the secondary-body-to-inertial DCM.
% The radial direction points from Didymos to Dimorphos. At synchronous
% alignment, +X_B=-e_r, +Y_B=-e_theta, and +Z_B=+h.
% The supplied orbit basis defines the inertial frame.

orbit_basis = [e_r0(:), e_th0(:), h0_hat(:)];
angle = theta + phi_B + pi;
rotation_body_to_orbit = [cos(angle), -sin(angle), 0; ...
    sin(angle),  cos(angle), 0; ...
    0,         0,        1];
rotation_body_to_inertial = orbit_basis * rotation_body_to_orbit;
end
