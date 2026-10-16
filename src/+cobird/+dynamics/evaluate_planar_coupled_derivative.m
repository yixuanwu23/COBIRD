function dydt = evaluate_planar_coupled_derivative(t, y, params, forces, F_Yark_planar, M_YORP_z)
%EVALUATE_PLANAR_COUPLED_DERIVATIVE Evaluate the coupled planar ODE.
% F_Yark_planar=[F_radial;F_tangential] is the in-plane recoil force;
% M_YORP_z is the secondary torque along the orbit normal.
% The recoil terms enter as
%   r̈   += F_radial     / m_Dim
%   θ̈   += F_tangential / (m_Dim * r)
%   φ̈_B += M_YORP_z    / I_{B,z}

% Unpack the state.
r      = y(1);
rdot   = y(2);
% theta = y(3);   % The inertial orbit angle is cyclic.
thdot  = y(4);
phi_B  = y(5);
phi_b_dot = y(6);

% Read precomputed parameters.
mu       = params.mu_system;
a2       = params.a2;
b2       = params.b2;
J2A_R2A  = params.J2A_R2A;
J2R2_Dim = params.J2R2_Dim;
C22R2    = params.C22R2;

mass_radius_inertia_ratio = 5 * params.mu_Didymos * r^2 / (mu * (a2 + b2));

cos_2phi = cos(2 * phi_B);
sin_2phi = sin(2 * phi_B);

% Apply force switches.
term_J2A  = 0;  if forces.J2_Didymos,    term_J2A  = J2A_R2A;  end
term_J2B  = 0;  if forces.J2_Dimorphos,  term_J2B  = J2R2_Dim; end
term_C22  = 0;  if forces.C22_Dimorphos, term_C22  = C22R2;    end

% Radial equation.
r_ddot = thdot^2 * r ...
    - (mu / r^2) * (1 + (3 / (2*r^2)) * ...
    (term_J2A + term_J2B + 6 * term_C22 * cos_2phi));

% Orbit-angle equation.
theta_ddot_drag = 0;
if forces.DART_drag
    et = params.et0 + t;
    if et >= params.t_imp
        dt_drag = et - params.t_imp;
        if dt_drag < params.drag_cutoff
            theta_ddot_drag = -params.drag_A * exp(-dt_drag / params.drag_tau);
        end
    end
end

th_ddot = -2 * rdot * thdot / r ...
    + (6 * mu * term_C22 * sin_2phi) / r^5 ...
    + theta_ddot_drag;

% Libration-angle equation.
phi_b_ddot = -(1 + mass_radius_inertia_ratio) * ...
    (6 * mu * term_C22 * sin_2phi) / r^5 ...
    + 2 * rdot * thdot / r;

% Naidu et al. (2024) Eq. (4): the fitted DART draglike term is kept
% as an empirical orbital acceleration only.  Real radiative force
% below remains torque-consistent and therefore does subtract its
% theta acceleration from the relative angle equation.

% Tidal torque.
if forces.tidal
    Gamma_B_mag = params.tidal.coeff_B / r^6 * params.tidal.k2_over_Q_secondary;
    delta_B = cobird.dynamics.compute_tidal_linearization_rate(Gamma_B_mag, params.tidal.I_Bz, ...
        params.tidal.rho_secondary, params.G_const);
    if Gamma_B_mag == 0
        sign_B = 0;
    elseif abs(phi_b_dot) > delta_B
        sign_B = sign(phi_b_dot);
    else
        sign_B = phi_b_dot / delta_B;
    end

    Gamma_B = -Gamma_B_mag * sign_B;
    Gamma_A = -params.tidal.coeff_A / r^6 * params.tidal.k2_over_Q_primary;

    theta_ddot_tidal = (-Gamma_A - Gamma_B) / (params.tidal.m_red * r^2);
    phi_b_ddot = phi_b_ddot + Gamma_B / params.tidal.I_Bz - theta_ddot_tidal;
    th_ddot  = th_ddot + theta_ddot_tidal;
end

% Yarkovsky force and YORP torque.
if forces.Yarkovsky
    m_Dim = params.mu_Dimorphos / params.G_const;
    I_Bz  = m_Dim * (a2 + b2) / 5;

    r_ddot   = r_ddot   + F_Yark_planar(1) / m_Dim;
    theta_ddot_yark = F_Yark_planar(2) / (m_Dim * r);
    th_ddot  = th_ddot + theta_ddot_yark;
    phi_b_ddot = phi_b_ddot + M_YORP_z / I_Bz - theta_ddot_yark;
end

% Assemble the state derivative.
dydt = [rdot; r_ddot; thdot; th_ddot; phi_b_dot; phi_b_ddot];

end
