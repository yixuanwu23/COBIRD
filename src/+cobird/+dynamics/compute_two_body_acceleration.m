function a = compute_two_body_acceleration(r, GM)
%COMPUTE_TWO_BODY_ACCELERATION Evaluate point-mass relative acceleration.
%   r is the secondary position relative to the primary.

a = -GM * r / norm(r)^3;
end
