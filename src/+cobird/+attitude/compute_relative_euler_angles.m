function [theta1, theta2, theta3] = compute_relative_euler_angles(r_J, v_J, q_B_to_J)
%COMPUTE_RELATIVE_EULER_ANGLES Dimorphos attitude relative to sync frame P.
% R_B_to_J = R_P_to_J * Rz(theta3) * Ry(theta2) * Rx(theta1), where
% R_P_to_J = [-e_r, -e_theta, h_hat]. Angles are returned in radians.
R_B_to_J = cobird.attitude.convert_quaternion_to_dcm(q_B_to_J);
[R_RTN_to_J, ~, ~, ~, ~] = cobird.attitude.build_rtn_frame(r_J, v_J);
R_P_to_J = R_RTN_to_J * diag([-1, -1, 1]);

% C maps P-frame components to B-frame components.
C_P_to_B = R_B_to_J' * R_P_to_J;
theta2 = asin(max(-1, min(1, -C_P_to_B(1,3))));
theta1 = atan2(C_P_to_B(2,3), C_P_to_B(3,3));
theta3 = atan2(C_P_to_B(1,2), C_P_to_B(1,1));
end
