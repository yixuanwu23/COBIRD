function q = convert_dcm_to_quaternion(R)
%CONVERT_DCM_TO_QUATERNION Convert a DCM to a scalar-first unit quaternion.
% Columns of R express source-frame axes in the target frame. The returned
% q=[q0;q1;q2;q3] follows the convention used by the attitude equations.

assert(isequal(size(R), [3, 3]) && all(isfinite(R), 'all'), ...
    'cobird:attitude:convert_dcm_to_quaternion:invalidDimensions', ...
    'The direction cosine matrix must be a finite 3-by-3 matrix.');
assert(norm(R' * R - eye(3), 'fro') < 1e-10 && det(R) > 0, ...
    'cobird:attitude:convert_dcm_to_quaternion:invalidRotation', ...
    'The input must be a right-handed orthogonal rotation matrix.');

candidates = [1 + R(1,1) + R(2,2) + R(3,3); ...
    1 + R(1,1) - R(2,2) - R(3,3); ...
    1 - R(1,1) + R(2,2) - R(3,3); ...
    1 - R(1,1) - R(2,2) + R(3,3)];
[~, idx] = max(candidates);

switch idx
    case 1
        q0 = 0.5 * sqrt(max(candidates(1), 0));
        q = [q0; (R(3,2)-R(2,3))/(4*q0); ...
            (R(1,3)-R(3,1))/(4*q0); ...
            (R(2,1)-R(1,2))/(4*q0)];
    case 2
        q1 = 0.5 * sqrt(max(candidates(2), 0));
        q = [(R(3,2)-R(2,3))/(4*q1); q1; ...
            (R(1,2)+R(2,1))/(4*q1); (R(1,3)+R(3,1))/(4*q1)];
    case 3
        q2 = 0.5 * sqrt(max(candidates(3), 0));
        q = [(R(1,3)-R(3,1))/(4*q2); ...
            (R(1,2)+R(2,1))/(4*q2); q2; (R(2,3)+R(3,2))/(4*q2)];
    otherwise
        q3 = 0.5 * sqrt(max(candidates(4), 0));
        q = [(R(2,1)-R(1,2))/(4*q3); ...
            (R(1,3)+R(3,1))/(4*q3); (R(2,3)+R(3,2))/(4*q3); q3];
end

q = q / norm(q);
if q(1) < 0
    q = -q;
end
end
