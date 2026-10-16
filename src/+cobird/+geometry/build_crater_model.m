function crater = build_crater_model(save_path, modelname, opts)
%BUILD_CRATER_MODEL Build the canonical hemispherical roughness model.
%
%   crater = cobird.geometry.build_crater_model(save_path, modelname, opts)
%
% The crater is a unit-radius hemispherical bowl whose opening is the
% polygonal aperture of the selected mesh.
%
% Fields:
%   vertices       [Nv x 3] local coordinates, z is parent normal direction
%   faces          [Nc x 3] triangulation, default Nc=100
%   normal         [Nc x 3] local emission normals into the crater void
%   area           [Nc x 1] unit-crater triangle areas
%   view           [Nc x Nc] M2 view factor from facet i to facet j
%   shadow_table   [Nc x N_bins] true if shadowed for a local sun bin

if nargin < 1 || isempty(save_path)
    save_path = pwd;
end
if nargin < 2 || isempty(modelname)
    modelname = 'hemisphere_crater_100facets';
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'N_az_mesh')
    opts.N_az_mesh = 20;
end
if ~isfield(opts, 'N_radial')
    opts.N_radial = 4;
end
if ~isfield(opts, 'N_az_shadow')
    opts.N_az_shadow = 360;
end
if ~isfield(opts, 'N_el_shadow')
    opts.N_el_shadow = 180;
end
if ~isfield(opts, 'ray_offset')
    opts.ray_offset = 1e-6;
end
if ~isfield(opts, 'opening_area_mode')
    opts.opening_area_mode = 'mesh_aperture';
end
opts.opening_area_mode = validatestring(opts.opening_area_mode, ...
    {'mesh_aperture'}, mfilename, 'opts.opening_area_mode');

validateattributes(opts.N_az_mesh, {'numeric'}, ...
    {'scalar', 'integer', '>=', 3}, mfilename, 'opts.N_az_mesh');
validateattributes(opts.N_radial, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, mfilename, 'opts.N_radial');
validateattributes(opts.N_az_shadow, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1}, mfilename, 'opts.N_az_shadow');
validateattributes(opts.N_el_shadow, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1}, mfilename, 'opts.N_el_shadow');
validateattributes(opts.ray_offset, {'numeric'}, ...
    {'scalar', 'positive', 'finite'}, mfilename, 'opts.ray_offset');

schema = 'COBIRD_HEMISPHERICAL_CRATER';
outfile = fullfile(save_path, [modelname, '.mat']);
if isfile(outfile)
    s = load(outfile, 'crater');
    crater = s.crater;
    if isfield(crater, 'schema') && strcmp(crater.schema, schema) && ...
            isfield(crater, 'opts') && crater.opts.N_az_mesh == opts.N_az_mesh && ...
            crater.opts.N_radial == opts.N_radial && ...
            crater.opts.N_az_shadow == opts.N_az_shadow && ...
            crater.opts.N_el_shadow == opts.N_el_shadow && ...
            crater.opts.ray_offset == opts.ray_offset && ...
            isfield(crater.opts, 'opening_area_mode') && ...
            strcmp(crater.opts.opening_area_mode, 'mesh_aperture')
        rim = crater.vertices(abs(crater.vertices(:, 3)) < 1e-10, 1:2);
        [~, order] = sort(atan2(rim(:, 2), rim(:, 1)));
        aperture = polyarea(rim(order, 1), rim(order, 2));
        assert(abs(crater.opening_area - aperture) < 1e-12 * aperture, ...
            'cobird:geometry:build_crater_model:aperture', ...
            'The cached opening area does not match the mesh aperture.');
        return
    end
    error('cobird:geometry:build_crater_model:incompatibleCache', ...
        ['The crater cache is incompatible with the requested geometry. ' ...
         'Rebuild it at a new path: %s'], outfile);
end

Naz = opts.N_az_mesh;
Nr  = opts.N_radial;
phi = (0:Naz-1)' * (2*pi/Naz);
theta = linspace(0, pi/2, Nr);  % bottom to rim

vertices = zeros(1 + (Nr-1)*Naz, 3);
vertices(1, :) = [0, 0, -1];
row = 2;
for r = 2:Nr
    th = theta(r);
    ring = [sin(th)*cos(phi), sin(th)*sin(phi), -cos(th)*ones(Naz, 1)];
    vertices(row:row+Naz-1, :) = ring;
    row = row + Naz;
end

faces = zeros(Naz + (Nr-2)*2*Naz, 3);
row = 1;
first_ring = 2:(Naz+1);
for a = 1:Naz
    b = mod(a, Naz) + 1;
    faces(row, :) = [1, first_ring(b), first_ring(a)];
    row = row + 1;
end
for r = 2:Nr-1
    ring0 = 2 + (r-2)*Naz;
    ring1 = 2 + (r-1)*Naz;
    for a = 1:Naz
        b = mod(a, Naz) + 1;
        v00 = ring0 + a - 1;
        v01 = ring0 + b - 1;
        v10 = ring1 + a - 1;
        v11 = ring1 + b - 1;
        faces(row, :) = [v00, v11, v10];
        row = row + 1;
        faces(row, :) = [v00, v01, v11];
        row = row + 1;
    end
end

[normal, ~, centroid] = local_geometry(vertices, faces);

% For an inner hemispherical surface, the emission normal points into
% the void. The triangle winding above is checked and flipped if needed.
radial_to_void = -centroid ./ vecnorm(centroid, 2, 2);
flip = dot(normal, radial_to_void, 2) < 0;
faces(flip, [2 3]) = faces(flip, [3 2]);
[normal, area, centroid] = local_geometry(vertices, faces);

Nc = size(faces, 1);
exchange_area = zeros(Nc, Nc);
view = zeros(Nc, Nc);
for i = 1:Nc-1
    for j = i+1:Nc
        displacement = centroid(j, :) - centroid(i, :);
        distance = norm(displacement);
        direction = displacement / distance;
        cosine_i = dot(normal(i, :), direction);
        cosine_j = dot(normal(j, :), -direction);
        if cosine_i <= 0 || cosine_j <= 0
            continue
        end
        term_i = atan(sqrt(pi * area(i)) * cosine_i / (2 * distance));
        term_j = atan(sqrt(pi * area(j)) * cosine_j / (2 * distance));
        exchange = 4 * sqrt(area(i) * area(j)) / pi^2 * term_i * term_j;
        exchange_area(i, j) = exchange;
        exchange_area(j, i) = exchange;
        view(i, j) = exchange / area(i);
        view(j, i) = exchange / area(j);
    end
end

row_sum = sum(view, 2);
assert(all(row_sum <= 1 + 1e-10), ...
    'cobird:geometry:build_crater_model:rowSumExceeded', ...
    'Crater M2 view-factor row sum exceeds one (max %.12g).', max(row_sum));
exchange_from_view = area .* view;
reciprocity_error = max(abs(exchange_from_view - exchange_from_view.'), [], 'all');
assert(reciprocity_error <= 1e-12 * max(area), ...
    'cobird:geometry:build_crater_model:reciprocityFailure', ...
    'Crater view factors violate reciprocity.');

shadow_table = local_shadow_table(vertices, faces, normal, centroid, opts);

rim = vertices(end-Naz+1:end, 1:2);
rim_next = rim([2:end, 1], :);
opening_area = 0.5 * abs(sum(rim(:, 1) .* rim_next(:, 2) ...
    - rim(:, 2) .* rim_next(:, 1)));

crater = struct();
crater.schema = schema;
crater.view_factor_method = 'M2_REZAC_ZHAO_2020';
crater.vertices = vertices;
crater.faces = faces;
crater.normal = normal;
crater.area = area;
crater.centroid = centroid;
crater.view = sparse(view);
crater.exchange_area = sparse(exchange_area);
crater.shadow_table = shadow_table;
crater.N_az_shadow = opts.N_az_shadow;
crater.N_el_shadow = opts.N_el_shadow;
crater.daz = 2*pi / opts.N_az_shadow;
crater.del = pi / opts.N_el_shadow;
crater.opening_area = opening_area;
crater.opening_area_mode = opts.opening_area_mode;
crater.opening_area_ratio_to_pi = opening_area / pi;
crater.sky_fraction = max(0, 1 - row_sum);
escape_area = area .* crater.sky_fraction;
crater.escape_area_ratio = sum(escape_area) / crater.opening_area;
crater.max_view_factor_row_sum = max(row_sum);
crater.reciprocity_error = reciprocity_error;
crater.opts = opts;

save(outfile, 'crater', '-v7.3');
end

function [normal, area, centroid] = local_geometry(vertices, faces)
e1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
e2 = vertices(faces(:, 3), :) - vertices(faces(:, 1), :);
cv = cross(e1, e2, 2);
norm_cv = vecnorm(cv, 2, 2);
normal = cv ./ norm_cv;
area = 0.5 * norm_cv;
centroid = (vertices(faces(:, 1), :) + vertices(faces(:, 2), :) + vertices(faces(:, 3), :)) / 3;
end

function shadow_table = local_shadow_table(vertices, faces, normal, centroid, opts)
if mod(opts.N_az_shadow, opts.N_az_mesh) == 0
    shadow_table = local_shadow_table_symmetric( ...
        vertices, faces, normal, centroid, opts);
else
    shadow_table = local_shadow_table_direct( ...
        vertices, faces, normal, centroid, opts);
end
end

function shadow_table = local_shadow_table_direct(vertices, faces, normal, centroid, opts)
Nc = size(faces, 1);
N_bins = opts.N_az_shadow * opts.N_el_shadow;
shadow_table = false(Nc, N_bins);

v1 = vertices(faces(:, 1), :);
v2 = vertices(faces(:, 2), :);
v3 = vertices(faces(:, 3), :);
for b = 1:N_bins
    az_idx = mod(b - 1, opts.N_az_shadow) + 1;
    el_idx = floor((b - 1) / opts.N_az_shadow) + 1;
    az = -pi + (az_idx - 0.5) * (2*pi/opts.N_az_shadow);
    el = -pi/2 + (el_idx - 0.5) * (pi/opts.N_el_shadow);
    ce = cos(el);
    sun = [ce*cos(az), ce*sin(az), sin(el)];

    lit = find(normal * sun' > 0);
    for kk = 1:numel(lit)
        i = lit(kk);
        ray_o = centroid(i, :) + normal(i, :) * opts.ray_offset;
        blockers = true(Nc, 1);
        blockers(i) = false;
        hit = cobird.geometry.intersect_triangle_ray(ray_o, sun, ...
            v1(blockers, :), v2(blockers, :), v3(blockers, :), ...
            'lineType', 'ray', 'planeType', 'two sided');
        shadow_table(i, b) = any(hit);
    end
end
end

function shadow_table = local_shadow_table_symmetric(vertices, faces, normal, centroid, opts)
% A regular crater mesh is rotationally symmetric. Only the azimuth
% bins in one mesh sector require ray tracing; the rest are rotations.
Nc = size(faces, 1);
Naz_mesh = opts.N_az_mesh;
Naz_shadow = opts.N_az_shadow;
Nel_shadow = opts.N_el_shadow;
bins_per_sector = Naz_shadow / Naz_mesh;
shadow_table = false(Nc, Naz_shadow * Nel_shadow);

[flat_origin, flat_v1, flat_v2, flat_v3, flat_ray_o] = ...
    local_flat_blockers(vertices, faces, normal, centroid, opts.ray_offset);
face_permutation = zeros(Nc, Naz_mesh);
for sector_shift = 0:Naz_mesh-1
    face_permutation(:, sector_shift + 1) = ...
        local_rotated_face_indices(Naz_mesh, opts.N_radial, sector_shift);
end

daz = 2*pi / Naz_shadow;
del = pi / Nel_shadow;
for el_idx = 1:Nel_shadow
    el = -pi/2 + (el_idx - 0.5) * del;
    if el <= 0
        continue
    end
    cos_el = cos(el);
    for base_az_idx = 1:bins_per_sector
        az = -pi + (base_az_idx - 0.5) * daz;
        sun = [cos_el*cos(az), cos_el*sin(az), sin(el)];
        lit = normal * sun' > 0;
        selected = lit(flat_origin);
        base_shadow = false(Nc, 1);
        if any(selected)
            hit = cobird.geometry.intersect_triangle_ray( ...
                flat_ray_o(selected, :), sun, ...
                flat_v1(selected, :), flat_v2(selected, :), flat_v3(selected, :), ...
                'lineType', 'ray', 'planeType', 'two sided');
            base_shadow = accumarray(flat_origin(selected), double(hit), ...
                [Nc, 1], @max, 0) > 0;
        end

        for sector_shift = 0:Naz_mesh-1
            target_az_idx = base_az_idx + sector_shift * bins_per_sector;
            bin = (el_idx - 1) * Naz_shadow + target_az_idx;
            perm = face_permutation(:, sector_shift + 1);
            shadow_table(perm, bin) = base_shadow;
        end
    end
end
end

function [flat_origin, flat_v1, flat_v2, flat_v3, flat_ray_o] = ...
    local_flat_blockers(vertices, faces, normal, centroid, ray_offset)
Nc = size(faces, 1);
pair_count = Nc * (Nc - 1);
flat_origin = zeros(pair_count, 1);
flat_blocker = zeros(pair_count, 1);
row = 1;
all_faces = (1:Nc)';
for i = 1:Nc
    blockers = all_faces(all_faces ~= i);
    rows = row:row + Nc - 2;
    flat_origin(rows) = i;
    flat_blocker(rows) = blockers;
    row = row + Nc - 1;
end
flat_v1 = vertices(faces(flat_blocker, 1), :);
flat_v2 = vertices(faces(flat_blocker, 2), :);
flat_v3 = vertices(faces(flat_blocker, 3), :);
ray_origin = centroid + normal .* ray_offset;
flat_ray_o = ray_origin(flat_origin, :);
end

function perm = local_rotated_face_indices(Naz, Nr, sector_shift)
Nc = Naz * (2 * Nr - 3);
perm = zeros(Nc, 1);
for az_idx = 1:Naz
    rotated_az = mod(az_idx - 1 + sector_shift, Naz) + 1;
    perm(az_idx) = rotated_az;
end
for band = 1:Nr-2
    base = Naz + (band - 1) * 2 * Naz;
    for az_idx = 1:Naz
        rotated_az = mod(az_idx - 1 + sector_shift, Naz) + 1;
        source = base + 2 * (az_idx - 1) + (1:2);
        target = base + 2 * (rotated_az - 1) + (1:2);
        perm(source) = target;
    end
end
end
