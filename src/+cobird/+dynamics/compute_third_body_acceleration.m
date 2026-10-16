function a = compute_third_body_acceleration(r, s, GM)
%COMPUTE_THIRD_BODY_ACCELERATION Return differential point-mass gravity.
% r and s are the target and perturber positions relative to the reference
% body [m]. GM is the perturber gravitational parameter [m^3/s^2].

% Target position relative to the perturber.
d = r - s;

% Subtract the acceleration of the reference body.
a = -GM * ( d/(norm(d)^3) + s/(norm(s)^3) );

end
