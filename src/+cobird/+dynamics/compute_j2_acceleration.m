function a_rel = compute_j2_acceleration(r, z_hat, mu_self, mu_other, J2R2)
%COMPUTE_J2_ACCELERATION Return the differential J2 acceleration.
% r and z_hat are expressed in ECLIPJ2000. J2R2 is J2*R_ref^2 [m^2].
% Newton's third law gives a_rel=(1+mu_other/mu_self)*a_J2.

rn = norm(r);
s  = dot(r, z_hat);
k  = 1.5 * mu_self * J2R2 / rn^5;

a_J2 = k * (r * (5*s^2/(rn*rn) - 1) - 2*s*z_hat);

a_rel = (1 + mu_other / mu_self) * a_J2;

end
