function w_primary_J2000 = compute_didymos_spin_vector(et, period_h)
%COMPUTE_DIDYMOS_SPIN_VECTOR Didymos PCK spin vector in ECLIPJ2000.
% et is SPICE ephemeris time [s]; period_h is the rotation period [h].
% The spin axis is +Z of IAU_DIDYMOS. The returned w_primary_J2000
% contains ECLIPJ2000 components [rad/s].

arguments
    et (1,1) double {mustBeFinite}
    period_h (1,1) double {mustBePositive, mustBeFinite}
end

R_Did_to_J2000 = cspice_pxform('IAU_DIDYMOS', 'ECLIPJ2000', et);
z_hat_J2000 = R_Did_to_J2000(:,3);
z_hat_J2000 = z_hat_J2000 / norm(z_hat_J2000);
w_primary_J2000 = z_hat_J2000 * (2*pi / (period_h * 3600));
end
