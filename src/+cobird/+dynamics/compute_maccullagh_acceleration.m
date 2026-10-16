function acceleration_inertial = compute_maccullagh_acceleration(r_hat_body, rn, I_hat, mu_system, rotation_body_to_inertial)
%COMPUTE_MACCULLAGH_ACCELERATION Return the non-point-mass MacCullagh term.
% r_hat_body and I_hat use the secondary body frame; I_hat is inertia per
% unit secondary mass [m^2]. The supplied rotation defines the inertial
% frame of the returned perturbing acceleration [m/s^2].
%
%   a_body = -(3*mu/(2*r^4)) * ...
%       (trace(Ihat)*rhat + 2*Ihat*rhat - 5*(rhat'*Ihat*rhat)*rhat)

tr_I = trace(I_hat);
i_times_r_hat = I_hat * r_hat_body;
r_i_r = r_hat_body' * i_times_r_hat;

acceleration_body = -(3 * mu_system / (2 * rn^4)) * ...
    (tr_I * r_hat_body - 5 * r_i_r * r_hat_body + 2 * i_times_r_hat);

acceleration_inertial = rotation_body_to_inertial * acceleration_body;

end
