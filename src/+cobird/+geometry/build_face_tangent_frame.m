function [x_local, y_local] = build_face_tangent_frame(vertices, faces, normal)
%BUILD_FACE_TANGENT_FRAME Build an orthonormal tangent basis on each face.

edge_1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
x_local = edge_1 ./ vecnorm(edge_1, 2, 2);
y_local = cross(normal, x_local, 2);
y_local = y_local ./ vecnorm(y_local, 2, 2);
end
