function system = make_synthetic_binary()
%MAKE_SYNTHETIC_BINARY Create a deterministic, data-free binary system.

system = struct();
system.primary = local_make_body([250, 220, 190], 2000);
system.secondary = local_make_body([100, 80, 65], 2000);
system.separation_m = 900;
end

function body = local_make_body(axes_m, density_kg_m3)
golden_ratio = (1 + sqrt(5)) / 2;
vertices = [ ...
    -1,  golden_ratio, 0;  1,  golden_ratio, 0; ...
    -1, -golden_ratio, 0;  1, -golden_ratio, 0; ...
     0, -1,  golden_ratio; 0, 1,  golden_ratio; ...
     0, -1, -golden_ratio; 0, 1, -golden_ratio; ...
     golden_ratio, 0, -1;  golden_ratio, 0, 1; ...
    -golden_ratio, 0, -1; -golden_ratio, 0, 1];
vertices = vertices ./ vecnorm(vertices, 2, 2);
vertices = vertices .* axes_m;

faces = [ ...
    1, 12, 6; 1, 6, 2; 1, 2, 8; 1, 8, 11; 1, 11, 12; ...
    2, 6, 10; 6, 12, 5; 12, 11, 3; 11, 8, 7; 8, 2, 9; ...
    4, 10, 5; 4, 5, 3; 4, 3, 7; 4, 7, 9; 4, 9, 10; ...
    5, 10, 6; 3, 5, 12; 7, 3, 11; 9, 7, 8; 10, 9, 2];

[normal, ~, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);
inward = sum(normal .* centroid, 2) < 0;
faces(inward, [2, 3]) = faces(inward, [3, 2]);
[normal, area, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);
[x_local, y_local] = cobird.geometry.build_face_tangent_frame(vertices, faces, normal);

face_count = size(faces, 1);
reabsorb = struct();
reabsorb.method = 'triangle_horizon';
reabsorb.face_occlusion_info = cell(face_count, 1);
reabsorb.candidate_faces = cell(face_count, 1);
reabsorb.has_candidates = false(face_count, 1);
no_occultation = true(face_count, 1);
occultation = cobird.geometry.build_occultation_precomputation( ...
    vertices, faces, normal, no_occultation, reabsorb);

shadow = struct();
shadow.N_az = 1;
shadow.N_el = 1;
shadow.daz = 2 * pi;
shadow.del = pi;
shadow.occ_rev_map = zeros(face_count, 1, 'uint32');
shadow.shadow_table = false(0, 1);

volume_m3 = 4 * pi / 3 * prod(axes_m);
body = struct();
body.axes_m = axes_m;
body.radius_m = max(axes_m);
body.density_kg_m3 = density_kg_m3;
body.mass_kg = density_kg_m3 * volume_m3;
body.vertices = vertices;
body.faces = faces;
body.normal = normal;
body.area = area;
body.centroid = centroid;
body.x_local = x_local;
body.y_local = y_local;
body.reabsorb = reabsorb;
body.occultation = occultation;
body.shadow = shadow;
end
