function [normal, area, centroid] = compute_face_geometry(vertices, faces)
%COMPUTE_FACE_GEOMETRY Compute unit normals, areas, and centroids of triangles.

edge_1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
edge_2 = vertices(faces(:, 3), :) - vertices(faces(:, 1), :);
cross_product = cross(edge_1, edge_2, 2);
double_area = vecnorm(cross_product, 2, 2);
normal = cross_product ./ double_area;
area = double_area / 2;
centroid = (vertices(faces(:, 1), :) + vertices(faces(:, 2), :) + ...
    vertices(faces(:, 3), :)) / 3;
end
