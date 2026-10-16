function R = convert_quaternion_to_dcm(q)
%CONVERT_QUATERNION_TO_DCM Convert a scalar-first quaternion to a DCM.
% For a body-to-J2000 quaternion q=[q0;q1;q2;q3], v_J2000=R*v_B.

q = q(:);
assert(numel(q) == 4 && all(isfinite(q)), ...
    'cobird:attitude:convert_quaternion_to_dcm:invalidDimensions', ...
    'The quaternion must be a finite four-element vector.');
q_norm = norm(q);
assert(q_norm > 0, 'cobird:attitude:convert_quaternion_to_dcm:zeroNorm', ...
    'The quaternion norm must be positive.');
q = q / q_norm;

q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
R = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3),     2*(q1*q3+q0*q2); ...
    2*(q1*q2+q0*q3),     q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1); ...
    2*(q1*q3-q0*q2),     2*(q2*q3+q0*q1),     q0^2-q1^2-q2^2+q3^2];
end
