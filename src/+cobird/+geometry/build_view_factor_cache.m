function build_view_factor_cache(vertices, faces, normal_all, area_all, centroid_all, ...
    output_dir, modelname, parameters)
%BUILD_VIEW_FACTOR_CACHE Build the sparse M2 radiative-exchange cache.

face_count = size(faces, 1);
if nargin < 8 || isempty(parameters)
    parameters = struct();
end
if ~isfield(parameters, 'F_threshold')
    parameters.F_threshold = 1e-4;
end
if ~isfield(parameters, 'ray_offset_m')
    parameters.ray_offset_m = 1e-6;
end

F_threshold = parameters.F_threshold;
ray_offset_m = parameters.ray_offset_m;
assert(isscalar(F_threshold) && isfinite(F_threshold) && F_threshold >= 0, ...
    'cobird:geometry:build_view_factor_cache:badThreshold', ...
    'F_threshold must be a finite nonnegative scalar.');
assert(isscalar(ray_offset_m) && isfinite(ray_offset_m) && ray_offset_m > 0, ...
    'cobird:geometry:build_view_factor_cache:badRayOffset', ...
    'ray_offset_m must be a finite positive scalar.');

fprintf('  [VF] %s: %d faces, M2, F_thr=%.1e\n', ...
    modelname, face_count, F_threshold);

pair_j_by_i = cell(face_count, 1);
exchange_area_by_i = cell(face_count, 1);

parfor i = 1:face_count
    normal_i = normal_all(i, :);
    centroid_i = centroid_all(i, :);
    area_i = area_all(i);
    all_candidates = cobird.geometry.find_potential_blockers( ...
        centroid_i, normal_i, vertices, faces);
    candidates = all_candidates(all_candidates > i);

    local_j = zeros(numel(candidates), 1);
    local_exchange_area = zeros(numel(candidates), 1);
    valid_count = 0;

    for candidate_index = 1:numel(candidates)
        j = candidates(candidate_index);
        normal_j = normal_all(j, :); %#ok<PFBNS>
        centroid_j = centroid_all(j, :); %#ok<PFBNS>
        area_j = area_all(j); %#ok<PFBNS>

        displacement = centroid_j - centroid_i;
        distance = norm(displacement);
        direction = displacement / distance;
        cosine_i = dot(normal_i, direction);
        cosine_j = dot(normal_j, -direction);
        if cosine_i <= 0 || cosine_j <= 0
            continue
        end

        exchange_area = local_m2_exchange_area( ...
            area_i, area_j, cosine_i, cosine_j, distance);
        if max(exchange_area / area_i, exchange_area / area_j) < F_threshold
            continue
        end

        blocker_index = all_candidates(all_candidates ~= i & all_candidates ~= j);
        if ~isempty(blocker_index)
            ray_origin = centroid_i + normal_i * ray_offset_m;
            ray_direction = centroid_j - ray_origin;
            blocked = cobird.geometry.intersect_triangle_ray( ...
                ray_origin, ray_direction, ...
                vertices(faces(blocker_index, 1), :), ...
                vertices(faces(blocker_index, 2), :), ...
                vertices(faces(blocker_index, 3), :), ...
                'lineType', 'segment', 'planeType', 'one sided');
            if any(blocked)
                continue
            end
        end

        valid_count = valid_count + 1;
        local_j(valid_count) = j;
        local_exchange_area(valid_count) = exchange_area;
    end

    pair_j_by_i{i} = local_j(1:valid_count);
    exchange_area_by_i{i} = local_exchange_area(1:valid_count);
end

counts = cellfun(@numel, pair_j_by_i);
pair_count = sum(counts);
pair_i = zeros(pair_count, 1);
pair_j = zeros(pair_count, 1);
exchange_area = zeros(pair_count, 1);

cursor = 1;
for i = 1:face_count
    count = counts(i);
    if count == 0
        continue
    end
    range = cursor:(cursor + count - 1);
    pair_i(range) = i;
    pair_j(range) = pair_j_by_i{i};
    exchange_area(range) = exchange_area_by_i{i};
    cursor = cursor + count;
end

exchange = sparse([pair_i; pair_j], [pair_j; pair_i], ...
    [exchange_area; exchange_area], face_count, face_count);
row_sum = full(sum(exchange, 2)) ./ area_all(:);
assert(all(row_sum <= 1 + 1e-10), ...
    'cobird:geometry:build_view_factor_cache:rowSumExceeded', ...
    'M2 view-factor row sum exceeds one (max %.12g).', max(row_sum));

VF = struct();
VF.schema = 'M2_VIEW_FACTOR_EXCHANGE';
VF.method = 'M2_REZAC_ZHAO_2020';
VF.pair_i = pair_i;
VF.pair_j = pair_j;
VF.exchange_area = exchange_area;
VF.N_faces = face_count;
VF.area = area_all(:);
VF.modelname = modelname;
VF.params = parameters;
VF.max_row_sum = max(row_sum);

output_path = fullfile(output_dir, [modelname, '_VF.mat']);
save(output_path, 'VF');
fprintf('  [VF] Saved %d pairs -> %s\n', pair_count, output_path);
end

function exchange_area = local_m2_exchange_area(area_i, area_j, cosine_i, cosine_j, distance)
term_i = atan(sqrt(pi * area_i) * cosine_i / (2 * distance));
term_j = atan(sqrt(pi * area_j) * cosine_j / (2 * distance));
exchange_area = 4 * sqrt(area_i * area_j) / pi^2 * term_i * term_j;
end
