function [tau_tidal_body, gamma_orb_J2000] = compute_tidal_response(r_J2000, v_J2000, ...
    w_body, R_B_to_J2000, params)
%COMPUTE_TIDAL_RESPONSE Return three-dimensional tidal torque and acceleration.
% Uses Cueva et al. (2024), Equations 24-26, with the three-dimensional
% direction from Meyer et al. (2023b), Equation 25, and the near-synchronous
% linearization from Jacobson and Scheeres (2011b), Appendix C.
%
% Inputs:
%   r_J2000, v_J2000  Secondary state relative to the primary [m, m/s]
%   w_body           Secondary angular velocity in its body frame [rad/s]
%   R_B_to_J2000      Body-to-inertial direction cosine matrix
%   params           Physical and tidal parameters
% All inertial vectors, including params.tidal.w_primary_J2000, must use
% the same frame. The output acceleration is expressed in that frame.
%
% Required parameter fields:
%   .tidal.k2_over_Q_secondary, .tidal.k2_over_Q_primary  Love number/Q
%   .tidal.w_primary_J2000  Primary inertial angular velocity [rad/s]
%   .Didymos_radius, .Dimorphos_radius  Tidal radii [m]
%   .Dim_a_ell, .Dim_b_ell  Secondary equatorial semiaxes [m]
%   .mu_Didymos, .mu_Dimorphos  Gravitational parameters [m^3/s^2]
%   .G_const  Gravitational constant [m^3/(kg s^2)]
% Didymos_radius sets the primary tidal torque and spherical inertia;
% R_Didymos is the separate J2 reference radius used by the orbit equations.
%
% Outputs:
%   tau_tidal_body   Specific Dimorphos torque in the body frame [m^2/s^2]
%   gamma_orb_J2000  Tidal orbital acceleration in the input frame [m/s^2]
%
% Torque magnitude (Murray and Dermott, 1999):
%    |Γ_i| = (3/2) G M_j² R_i⁵ / r⁶ × (k₂/Q)_i
%
% Three-dimensional direction:
%    Γ̂ = −φ̇_⊥ / |φ̇_⊥|
%    φ̇_⊥ = (ω_body − ω_orb) − [(ω_body − ω_orb)·r̂] r̂
%
% Near-synchronous linearization for |φ̇_⊥| < ε:
%    Γ̂ ≈ −φ̇_⊥ / ε
%
% Orbital acceleration from angular-momentum conservation:
%    Γ_orb = −(Γ_A + Γ_B)
%    γ_orb = (Γ_orb × r) / (m_red × r²)

G = params.G_const;

% Relative geometry and orbital angular velocity.
rn          = norm(r_J2000);
r_hat_J2000 = r_J2000 / rn;

% w_orb=(r x v)/r^2 points along the orbit normal.
w_orb_J2000 = cross(r_J2000, v_J2000) / rn^2;

% Component and reduced masses.
M_A   = params.mu_Didymos   / G;
M_B   = params.mu_Dimorphos / G;
m_red = M_A * M_B / (M_A + M_B);

R_A = params.Didymos_radius;
R_B = params.Dimorphos_radius;

% Tidal torque magnitudes [N m].
Gmag_B = 1.5 * G * M_A^2 * R_B^5 / rn^6 * params.tidal.k2_over_Q_secondary;
Gmag_A = 1.5 * G * M_B^2 * R_A^5 / rn^6 * params.tidal.k2_over_Q_primary;

% Dimorphos torque direction in the body frame.
r_hat_body = R_B_to_J2000' * r_hat_J2000;

% Express orbital angular velocity in the body frame.
w_orb_body = R_B_to_J2000' * w_orb_J2000;

% Remove the radial component of the relative spin.
phi_dot_B      = w_body - w_orb_body;
phi_dot_B_perp = phi_dot_B - dot(phi_dot_B, r_hat_body) * r_hat_body;
mag_B          = norm(phi_dot_B_perp);

I_Bz = M_B * (params.Dim_a_ell^2 + params.Dim_b_ell^2) / 5;
rho_B = M_B / ((4*pi/3) * R_B^3);
delta_B = cobird.dynamics.compute_tidal_linearization_rate(Gmag_B, I_Bz, rho_B, G);

if Gmag_B == 0
    Ghat_body = zeros(3,1);
elseif mag_B > delta_B
    Ghat_body = -phi_dot_B_perp / mag_B;
else
    Ghat_body = -phi_dot_B_perp / delta_B;
end

% Specific torque enters the attitude equation directly.
tau_tidal_body = (Gmag_B / M_B) * Ghat_body;

% Primary torque in the common inertial frame, used for orbital acceleration.
phi_dot_A      = params.tidal.w_primary_J2000 - w_orb_J2000;
phi_dot_A_perp = phi_dot_A - dot(phi_dot_A, r_hat_J2000) * r_hat_J2000;
mag_A          = norm(phi_dot_A_perp);

I_Az = (2/5) * M_A * R_A^2;
rho_A = M_A / ((4*pi/3) * R_A^3);
delta_A = cobird.dynamics.compute_tidal_linearization_rate(Gmag_A, I_Az, rho_A, G);

if Gmag_A == 0
    Ghat_A = zeros(3,1);
elseif mag_A > delta_A
    Ghat_A = -phi_dot_A_perp / mag_A;
else
    Ghat_A = -phi_dot_A_perp / delta_A;
end

Gamma_A_J2000 = Gmag_A * Ghat_A;   % Primary tidal torque [N m]

% Rotate the secondary torque to the common inertial frame.
Gamma_B_J2000 = Gmag_B * (R_B_to_J2000 * Ghat_body);

% Orbital acceleration (Cueva et al., 2024, Equation 26).
%   Γ_orb = −(Γ_A + Γ_B)
%   γ_orb = (Γ_orb × r) / (m_red × r²)
Gamma_orb       = -(Gamma_A_J2000 + Gamma_B_J2000);
gamma_orb_J2000 = cross(Gamma_orb, r_J2000) / (m_red * rn^2);

end
