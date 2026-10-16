function dydt = evaluate_3d_coupled_derivative( ...
    t, y, params, forces, F_Dim_body, F_Did_J2000, M_YORP_body, z_hat_Did_J2000, ...
    spice_Sun, spice_EMb, spice_Ven, spice_Jup)
%EVALUATE_3D_COUPLED_DERIVATIVE Evaluate the coupled thermal-orbit-attitude ODE.
% Includes Yarkovsky force and YORP torque when forces.Yarkovsky is true.
% The 13-state layout contains relative inertial position and velocity
% [m, m/s], a scalar-first body-to-inertial quaternion, and body-frame
% angular velocity [rad/s]. All inertial inputs must use the same frame.
% F_Dim_body and M_YORP_body use the secondary body frame [N, N m].
% F_Did_J2000 uses the inertial frame [N]; z_hat_Did_J2000 is the primary
% spin-axis unit vector in that frame. Enabled third-body positions are
% finite inertial column vectors relative to the primary [m].
% forces.secondary_gravity (default true) controls the secondary quadrupole
% and its gravity-gradient torque. The rigid-body inertia params.I is unchanged.

secondary_gravity = ~isfield(forces, 'secondary_gravity') || forces.secondary_gravity;

% Unpack the state.
r     = y(1:3);
v     = y(4:6);
q_raw = y(7:10);
w     = y(11:13);

% Normalize the quaternion and form the DCM.
q_norm_sq = dot(q_raw, q_raw);
q  = q_raw / sqrt(q_norm_sq);
b0 = q(1); b1 = q(2); b2 = q(3); b3 = q(4);

R_B_to_J2000 = cobird.attitude.convert_quaternion_to_dcm(q);

% Attitude dynamics.
rn = norm(r);

r_hat_body = R_B_to_J2000' * (r / rn);

tau_grav = zeros(3,1);
if secondary_gravity
    c_hat_body = -r_hat_body;
    tau_grav = (3 * params.mu_Didymos / rn^3) * ...
        cross(c_hat_body, params.I * c_hat_body);
end

tau_tidal       = zeros(3,1);
gamma_orb_tidal = zeros(3,1);
if forces.tidal
    [tau_tidal, gamma_orb_tidal] = cobird.dynamics.compute_tidal_response(r, v, w, R_B_to_J2000, params);
end

% Euler equation with optional YORP torque.
if forces.Yarkovsky
    m_Dim = params.mu_Dimorphos / params.G_const;
    tau_total = tau_grav + tau_tidal + M_YORP_body / m_Dim;
else
    tau_total = tau_grav + tau_tidal;
end
dw = params.I \ (tau_total - cross(w, params.I * w));

% Quaternion kinematics.
Xi_q = [-b1, -b2, -b3;
    b0, -b3,  b2;
    b3,  b0, -b1;
    -b2,  b1,  b0];
dq = 0.5 * Xi_q * w;

% Orbital dynamics.
a = zeros(3,1);

a = a + gamma_orb_tidal;
a = a + cobird.dynamics.compute_two_body_acceleration(r, params.mu_system);
if secondary_gravity
    a = a + cobird.dynamics.compute_maccullagh_acceleration(r_hat_body, rn, params.I, params.mu_system, R_B_to_J2000);
end

if forces.J2_Didymos
    a = a + cobird.dynamics.compute_j2_acceleration(r, z_hat_Did_J2000, params.mu_Didymos, params.mu_Dimorphos, ...
        params.J2_Didymos * params.R_Didymos^2);
end

% Third-body positions are sampled once per thermal step by the caller.
if forces.Sun
    a = a + cobird.dynamics.compute_third_body_acceleration(r, spice_Sun, params.GM_Sun);
end
if forces.EMB
    a = a + cobird.dynamics.compute_third_body_acceleration(r, spice_EMb, 403503.235507e9);
end
if forces.Venus
    a = a + cobird.dynamics.compute_third_body_acceleration(r, spice_Ven, 324858.592000e9);
end
if forces.Jupiter
    a = a + cobird.dynamics.compute_third_body_acceleration(r, spice_Jup, 126712764.100000e9);
end

et = params.et0 + t;

if forces.DART_drag
    a = a + cobird.dynamics.compute_dart_drag_acceleration(r, v, et, params);
end

% Rotate the secondary recoil force to the common inertial frame.
if forces.Yarkovsky
    m_Did = params.mu_Didymos / params.G_const;
    % Split-step contract over one thermal window:
    % - Dimorphos force is frozen in B and follows the integrated B(t).
    % - Didymos force is frozen in the inertial frame.
    a = a + (R_B_to_J2000 * F_Dim_body) / m_Dim ...
        - F_Did_J2000 / m_Did;
end

% Assemble the state derivative.
dydt = [v; a; dq; dw];

end
