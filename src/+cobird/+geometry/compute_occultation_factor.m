function f = compute_occultation_factor(r_didymos, r_Sun, R_sun_radius, R_did_radius)
%COMPUTE_OCCULTATION_FACTOR Return the visible fraction of the solar disk.
validateattributes(r_didymos, {'numeric'}, {'real','finite','vector','numel',3});
validateattributes(r_Sun, {'numeric'}, {'real','finite','vector','numel',3});
validateattributes(R_sun_radius, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(R_did_radius, {'numeric'}, {'real','finite','scalar','positive'});

d_did = norm(r_didymos);
d_sun = norm(r_Sun);
if d_did <= R_did_radius || d_sun <= R_sun_radius
    error('cobird:geometry:compute_occultation_factor:observerInsideBody', ...
        'The observer must lie outside both the Didymos and solar spheres.');
end

theta_did = asin(min(1, R_did_radius / d_did));
theta_sun = asin(min(1, R_sun_radius / d_sun));
cos_theta_sep = dot(r_didymos, r_Sun) / (d_did * d_sun);
theta_sep = acos(max(-1, min(1, cos_theta_sep)));

if theta_sep >= theta_did + theta_sun
    overlap_area = 0;
elseif theta_sep <= abs(theta_did - theta_sun)
    overlap_area = pi * min(theta_did, theta_sun)^2;
else
    arg_did = (theta_sep^2 + theta_did^2 - theta_sun^2) / ...
        (2 * theta_sep * theta_did);
    arg_sun = (theta_sep^2 + theta_sun^2 - theta_did^2) / ...
        (2 * theta_sep * theta_sun);
    radicand = (-theta_sep + theta_did + theta_sun) * ...
        (theta_sep + theta_did - theta_sun) * ...
        (theta_sep - theta_did + theta_sun) * ...
        (theta_sep + theta_did + theta_sun);
    overlap_area = theta_did^2 * acos(max(-1, min(1, arg_did))) + ...
        theta_sun^2 * acos(max(-1, min(1, arg_sun))) - ...
        0.5 * sqrt(max(0, radicand));
end

f = max(0, min(1, 1 - overlap_area / (pi * theta_sun^2)));
end
