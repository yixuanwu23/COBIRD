function entry = prepare_model_cache(model_root, model_name, mode, expected)
%PREPARE_MODEL_CACHE Check, complete, or rebuild one model cache set.
% model_name is an OBJ path stem relative to model_root. Vertex coordinates
% in this file must be in metres. For source meshes in other units, convert
% the loaded vertices to SI and pass those arrays to the lower-level cache
% builders; retain the original source mesh and its metadata unchanged.
% mode='check' validates existing products, 'ensure' builds missing ones,
% and 'rebuild' regenerates the cache set.

arguments
    model_root (1,:) char
    model_name (1,:) char
    mode (1,:) char {mustBeMember(mode, {'check', 'ensure', 'rebuild'})} = 'check'
    expected struct = struct('N_az', 360, 'N_el', 180, ...
        'F_threshold', 1e-4, 'ray_offset_m', 1e-6)
end

required_settings = {'N_az', 'N_el', 'F_threshold', 'ray_offset_m'};
assert(all(isfield(expected, required_settings)), ...
    'cobird:geometry:prepare_model_cache:missingSettings', ...
    'expected must contain N_az, N_el, F_threshold, and ray_offset_m.');

obj_path = fullfile(model_root, [model_name, '.obj']);
assert(isfile(obj_path), 'cobird:geometry:prepare_model_cache:missingObj', ...
    'Model OBJ file not found: %s', obj_path);
[model_dir, stem] = fileparts(obj_path);
paths = struct( ...
    'reabsorb', fullfile(model_dir, [stem, '_reabsorb_info.mat']), ...
    'shadow', fullfile(model_dir, [stem, '_shadow_horizon.mat']), ...
    'view_factor', fullfile(model_dir, [stem, '_VF.mat']));

manifest_entry = local_manifest_entry(model_root, model_name);
if ~isempty(manifest_entry) && ~strcmp(mode, 'rebuild')
    obj_hash = cobird.io.compute_file_sha256(obj_path);
    assert(strcmpi(obj_hash, manifest_entry.obj_sha256), ...
        'cobird:geometry:prepare_model_cache:objChanged', ...
        ['The OBJ file differs from the manifest, so caches will not be overwritten: %s. ' ...
         'After confirming the model change, rerun explicitly in rebuild mode.'], obj_path);
    local_reject_corrupt_existing(paths, manifest_entry);
end

missing = local_missing_products(paths);
if strcmp(mode, 'ensure') && isempty(manifest_entry) && numel(missing) < 3
    error('cobird:geometry:prepare_model_cache:untrackedProducts', ...
        ['The manifest has no model record, but caches already exist and cannot be verified: %s. ' ...
         'Inspect them and rerun explicitly in rebuild mode.'], model_name);
end
if strcmp(mode, 'check')
    assert(~isempty(manifest_entry), 'cobird:geometry:prepare_model_cache:missingManifestEntry', ...
        'The manifest has no record for model %s.', model_name);
    assert(isempty(missing), 'cobird:geometry:prepare_model_cache:missingProducts', ...
        'Model caches are missing: %s', strjoin(missing, ', '));
    entry = cobird.geometry.validate_manifest(model_root, model_name, expected);
    return
end

if strcmp(mode, 'ensure') && ~isempty(manifest_entry) && isempty(missing)
    entry = cobird.geometry.validate_manifest(model_root, model_name, expected);
    return
end

if strcmp(mode, 'rebuild')
    products_to_build = {'reabsorb', 'shadow', 'view_factor'};
elseif isempty(manifest_entry)
    products_to_build = {'reabsorb', 'shadow', 'view_factor'};
else
    products_to_build = missing;
end

[vertices, faces] = cobird.geometry.load_obj_mesh(obj_path);
[normal, area, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);
if ismember('reabsorb', products_to_build)
    fprintf('Building reabsorption cache for %s.\n', model_name);
    cobird.geometry.build_reabsorption_cache(vertices, faces, model_dir, stem);
end
if ismember('shadow', products_to_build)
    assert(isfile(paths.reabsorb), ...
        'cobird:geometry:prepare_model_cache:missingReabsorb', ...
        'A reabsorption cache is required before building the shadow cache.');
    data = load(paths.reabsorb, 'reabsorb_info');
    assert(isfield(data, 'reabsorb_info'), ...
        'cobird:geometry:prepare_model_cache:badReabsorb', ...
        'The reabsorption cache is incomplete.');
    fprintf('Building shadow cache for %s.\n', model_name);
    cobird.geometry.build_shadow_horizon_cache(vertices, faces, normal, centroid, ...
        data.reabsorb_info, expected.N_az, expected.N_el, model_dir, stem);
end
if ismember('view_factor', products_to_build)
    fprintf('Building view-factor cache for %s.\n', model_name);
    vf_parameters = struct('F_threshold', expected.F_threshold, ...
        'ray_offset_m', expected.ray_offset_m);
    cobird.geometry.build_view_factor_cache(vertices, faces, normal, area, centroid, ...
        model_dir, stem, vf_parameters);
end

local_validate_generated(paths, size(faces, 1), expected);
cobird.geometry.update_manifest(model_root, model_name);
entry = cobird.geometry.validate_manifest(model_root, model_name, expected);
end

function entry = local_manifest_entry(model_root, model_name)
manifest_path = fullfile(model_root, 'processed', 'metadata', 'geometry_manifest.json');
entry = [];
if ~isfile(manifest_path)
    return
end
manifest = jsondecode(fileread(manifest_path));
if ~isfield(manifest, 'models')
    return
end
normalized = strrep(model_name, '\', '/');
matches = strcmp({manifest.models.modelname}, normalized);
assert(nnz(matches) <= 1, 'cobird:geometry:prepare_model_cache:duplicateManifestEntry', ...
    'The manifest contains duplicate records for model %s.', model_name);
if any(matches)
    entry = manifest.models(matches);
end
end

function local_reject_corrupt_existing(paths, entry)
names = {'reabsorb', 'shadow', 'view_factor'};
hash_fields = {'reabsorb_sha256', 'shadow_sha256', 'view_factor_sha256'};
for k = 1:numel(names)
    path = paths.(names{k});
    if isfile(path)
        actual = cobird.io.compute_file_sha256(path);
        assert(strcmpi(actual, entry.(hash_fields{k})), ...
            'cobird:geometry:prepare_model_cache:cacheChanged', ...
            ['A cache exists but differs from the manifest and will not be overwritten: %s. ' ...
             'After confirming the change, rerun explicitly in rebuild mode.'], path);
    end
end
end

function missing = local_missing_products(paths)
names = {'reabsorb', 'shadow', 'view_factor'};
missing = names(~cellfun(@(name) isfile(paths.(name)), names));
end

function local_validate_generated(paths, face_count, expected)
assert(all(structfun(@isfile, paths)), ...
    'cobird:geometry:prepare_model_cache:generationIncomplete', ...
    'Cache generation did not produce all required files.');
reabsorb = load(paths.reabsorb, 'reabsorb_info');
assert(isfield(reabsorb, 'reabsorb_info') && ...
    numel(reabsorb.reabsorb_info.face_occlusion_info) == face_count && ...
    numel(reabsorb.reabsorb_info.candidate_faces) == face_count, ...
    'cobird:geometry:prepare_model_cache:badReabsorb', ...
    'The reabsorption cache has invalid dimensions.');
shadow = load(paths.shadow, 'schema', 'shadow_table', 'occ_rev_map', 'N_az', 'N_el');
assert(isfield(shadow, 'schema') && strcmp(shadow.schema, 'COBIRD_SHADOW_HORIZON') && ...
    shadow.N_az == expected.N_az && shadow.N_el == expected.N_el && ...
    numel(shadow.occ_rev_map) == face_count && ...
    size(shadow.shadow_table, 2) == expected.N_az * expected.N_el, ...
    'cobird:geometry:prepare_model_cache:badShadow', ...
    'The shadow cache has invalid dimensions or resolution.');
view_data = load(paths.view_factor, 'VF');
assert(isfield(view_data, 'VF') && ...
    isfield(view_data.VF, 'schema') && ...
    strcmp(view_data.VF.schema, 'M2_VIEW_FACTOR_EXCHANGE') && ...
    strcmp(view_data.VF.method, 'M2_REZAC_ZHAO_2020') && ...
    view_data.VF.N_faces == face_count && ...
    numel(view_data.VF.area) == face_count && ...
    all(view_data.VF.exchange_area >= 0) && ...
    abs(view_data.VF.params.F_threshold - expected.F_threshold) <= eps && ...
    abs(view_data.VF.params.ray_offset_m - expected.ray_offset_m) <= eps, ...
    'cobird:geometry:prepare_model_cache:badViewFactor', ...
    'The view-factor cache is invalid.');
end
