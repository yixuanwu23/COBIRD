function candidate_faces_index = find_potential_blockers(face_centroid, normal, vertices, faces)
%FIND_POTENTIAL_BLOCKERS Find faces above the tangent plane of a source facet.

% Vectorized: compute all dot products at once with matrix ops (no loop)

% Displacement from face_centroid to every vertex [V×3]
verts_shifted = vertices - face_centroid;          % broadcast [1×3]

% Normalize each row
norms_v = vecnorm(verts_shifted, 2, 2);            % [V×1]
norms_v(norms_v < eps) = 1;                        % guard against zero
verts_norm = verts_shifted ./ norms_v;             % [V×3]

% Dot product of every vertex direction with the face normal [V×1]
all_dots = verts_norm * normal';                   % normal is [1×3]

% For each face, check if ANY vertex has dot > threshold
dot1 = all_dots(faces(:,1));                       % [N×1]
dot2 = all_dots(faces(:,2));                       % [N×1]
dot3 = all_dots(faces(:,3));                       % [N×1]

candidate_faces_index = find(dot1 > 1e-8 | dot2 > 1e-8 | dot3 > 1e-8);
end
