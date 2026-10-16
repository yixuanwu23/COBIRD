function result = build_shadow_horizon_cache(vertices, faces, normal, centroid, ...
    reabsorb_info, N_az, N_el, save_path, model_name)
%BUILD_SHADOW_HORIZON_CACHE Build a body-frame self-shadow lookup table.
% For each discretized direction, mark facets blocked by local topography.
% The [N_occ x N_bins] table supports constant-time solar and mutual-
% radiation shadow queries without per-step ray tracing.
%
%   result = cobird.geometry.build_shadow_horizon_cache(vertices, faces, normal, centroid, ...
%       reabsorb_info, N_az, N_el, save_path, model_name)
%
% Inputs:
%   vertices       [Nv x 3] Mesh vertices in the body frame [m]
%   faces          [Nf x 3] Triangle vertex indices
%   normal         [Nf x 3] Unit facet normals
%   centroid       [Nf x 3] Facet centroids [m]
%   reabsorb_info  struct   Output from cobird.geometry.build_reabsorption_cache
%   N_az, N_el     scalar   Azimuth and elevation bin counts
%   save_path      string   Optional output directory
%   model_name     string   Model stem
%
% Output:
%   result  struct with fields:
%     shadow_table   [N_occ x N_bins] logical, true when blocked
%     occ_faces      [N_occ x 1] global facet indices
%     N_az, N_el     scalar bin counts
%     daz, del       scalar bin widths [rad]

assert(size(faces, 2) == 3, ...
    'cobird:geometry:build_shadow_horizon_cache:invalidFaces', ...
    'faces must be an N-by-3 array.');
if nargin < 6 || isempty(N_az)
    N_az = 360;
end
if nargin < 7 || isempty(N_el)
    N_el = 180;
end

EPS_RAY = 1e-2;   % Normal offset of each ray origin [m]

% Identify facets with potential blockers.
no_occ    = cellfun(@isempty, reabsorb_info.face_occlusion_info);
occ_faces = find(~no_occ);
N_occ     = numel(occ_faces);

if N_occ == 0
    result = struct('schema', 'COBIRD_SHADOW_HORIZON', 'ray_offset_m', EPS_RAY, ...
        'shadow_table', false(0, N_az * N_el), 'occ_faces', zeros(0, 1), ...
        'occ_rev_map', zeros(size(normal, 1), 1, 'uint32'), ...
        'N_az', N_az, 'N_el', N_el, 'daz', 2*pi/N_az, 'del', pi/N_el);
    if nargin >= 9 && ~isempty(save_path)
        outfile = fullfile(save_path, [model_name, '_shadow_horizon.mat']);
        save(outfile, '-struct', 'result', '-v7.3');
    end
    fprintf('No occluded faces; the shadow table is empty.\n');
    return
end

% Flatten candidate triangles for vectorized ray tracing.
M_arr  = cellfun(@(x) size(x, 1), reabsorb_info.candidate_faces(occ_faces));
N_flat = sum(M_arr);

flat_v1    = zeros(N_flat, 3);
flat_v2    = zeros(N_flat, 3);
flat_v3    = zeros(N_flat, 3);
flat_ray_o = zeros(N_flat, 3);
flat_fid   = zeros(N_flat, 1);   % 1-based index into occ_faces

row = 1;
for k = 1:N_occ
    fi  = occ_faces(k);
    cf  = reabsorb_info.candidate_faces{fi};
    M_k = size(cf, 1);
    rng = row : row + M_k - 1;
    flat_v1(rng, :)    = vertices(cf(:, 1), :);
    flat_v2(rng, :)    = vertices(cf(:, 2), :);
    flat_v3(rng, :)    = vertices(cf(:, 3), :);
    flat_ray_o(rng, :) = repmat(centroid(fi, :) + normal(fi, :) * EPS_RAY, M_k, 1);
    flat_fid(rng)     = k;
    row = row + M_k;
end

% Direction-bin parameters.
daz    = 2 * pi / N_az;
del    = pi / N_el;
N_bins = N_az * N_el;

shadow_table = false(N_occ, N_bins);

% Normals of potentially occluded facets.
n_occ = normal(occ_faces, :);   % [N_occ x 3]

fprintf('Shadow horizon: %d occ faces × %d bins (%d az × %d el), N_flat=%d\n', ...
    N_occ, N_bins, N_az, N_el, N_flat);
t0 = tic;

parfor b = 1:N_bins
    az_idx = mod(b - 1, N_az) + 1;
    el_idx = floor((b - 1) / N_az) + 1;
    az = -pi + (az_idx - 0.5) * daz;
    el = -pi/2 + (el_idx - 0.5) * del;

    % Body-frame unit direction at the bin center.
    cos_el = cos(el);
    d_sun  = [cos_el * cos(az), cos_el * sin(az), sin(el)];

    shadow_column = false(N_occ, 1);
    dot_vals = n_occ * d_sun';          % [N_occ x 1]
    lit_idx  = find(dot_vals > 0);
    if ~isempty(lit_idx)
        lit_set = false(N_occ, 1);
        lit_set(lit_idx) = true;
        selected = lit_set(flat_fid);
        if any(selected)
            remap = zeros(N_occ, 1);
            remap(lit_idx) = 1:numel(lit_idx);
            local_ids = remap(flat_fid(selected));
            hit = cobird.geometry.intersect_triangle_ray( ...
                flat_ray_o(selected, :), d_sun, ... %#ok<PFBNS>
                flat_v1(selected, :), flat_v2(selected, :), flat_v3(selected, :)); %#ok<PFBNS>
            any_hit = accumarray(local_ids, double(hit), ...
                [numel(lit_idx), 1]) > 0;
            shadow_column(lit_idx(any_hit)) = true;
        end
    end
    shadow_table(:, b) = shadow_column;
end

elapsed = toc(t0);
fprintf('Shadow horizon done: %.1f s\n', elapsed);
fprintf('Table: %d × %d, %.2f MB, density %.1f%%\n', ...
    N_occ, N_bins, N_occ * N_bins / 1e6, ...
    100 * nnz(shadow_table) / numel(shadow_table));

% Package the lookup table and metadata.
result = struct();
result.schema       = 'COBIRD_SHADOW_HORIZON';
result.ray_offset_m = EPS_RAY;
result.shadow_table = shadow_table;
result.occ_faces    = occ_faces;
result.N_az         = N_az;
result.N_el         = N_el;
result.daz          = daz;
result.del          = del;

% Map global facet indices to table rows; zero denotes no blockers.
N_total = size(normal, 1);
occ_rev_map = zeros(N_total, 1, 'uint32');
occ_rev_map(occ_faces) = uint32(1:N_occ);
result.occ_rev_map = occ_rev_map;

% Save when an output location is supplied.
if nargin >= 8 && ~isempty(save_path) && nargin >= 9
    outfile = fullfile(save_path, [model_name, '_shadow_horizon.mat']);
    save(outfile, '-struct', 'result', '-v7.3');
    fprintf('Saved: %s\n', outfile);
end

end
