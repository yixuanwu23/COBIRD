function build_reabsorption_cache(vertices, faces, output_dir, model_name)
%BUILD_REABSORPTION_CACHE Precompute terrain horizons for radiative recoil.
face_occlusion_info = cell(size(faces, 1), 1);
candidate_faces = cell(size(faces, 1), 1);

parfor i = 1:size(faces, 1)
    edge = vertices(faces(i, 2), :) - vertices(faces(i, 1), :); %#ok<PFBNS>
    other_edge = vertices(faces(i, 3), :) - vertices(faces(i, 1), :);
    normal = cross(edge, other_edge);
    normal = normal / norm(normal);
    centroid = mean(vertices(faces(i, :), :), 1);
    x_local = edge / norm(edge);
    y_local = cross(normal, x_local);
    frame = [x_local; y_local; normal]';

    indices = cobird.geometry.find_potential_blockers(centroid, normal, vertices, faces);
    if isempty(indices)
        continue
    end
    candidates = faces(indices, :);
    v1 = (vertices(candidates(:, 1), :) - centroid) * frame;
    v2 = (vertices(candidates(:, 2), :) - centroid) * frame;
    v3 = (vertices(candidates(:, 3), :) - centroid) * frame;

    % Keep the same near-horizon cutoff used for shadow-cache candidates.
    elevation = [atan2(v1(:, 3), hypot(v1(:, 1), v1(:, 2))), ...
        atan2(v2(:, 3), hypot(v2(:, 1), v2(:, 2))), ...
        atan2(v3(:, 3), hypot(v3(:, 1), v3(:, 2)))];
    keep = any(elevation >= deg2rad(0.01), 2);
    if any(keep)
        face_occlusion_info{i} = cobird.geometry.integrate_terrain_horizon( ...
            v1(keep, :), v2(keep, :), v3(keep, :));
        candidate_faces{i} = candidates(keep, :);
    end
end

reabsorb_info = struct('method', 'triangle_horizon', ...
    'face_occlusion_info', {face_occlusion_info}, ...
    'candidate_faces', {candidate_faces}, ...
    'has_candidates', ~cellfun(@isempty, candidate_faces));
save(fullfile(output_dir, [model_name, '_reabsorb_info.mat']), 'reabsorb_info');
end
