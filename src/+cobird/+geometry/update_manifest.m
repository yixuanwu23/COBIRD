function manifest = update_manifest(model_root, model_name)
%UPDATE_MANIFEST Record one newly generated model cache set.
% Other model records are preserved, including their original file hashes.

arguments
    model_root (1,:) char
    model_name (1,:) char
end

manifest_path = fullfile(model_root, 'processed', 'metadata', 'geometry_manifest.json');
if isfile(manifest_path)
    manifest = jsondecode(fileread(manifest_path));
    assert(isfield(manifest, 'schema') && strcmp(manifest.schema, 'COBIRD_GEOMETRY_MANIFEST') && ...
        isfield(manifest, 'models'), 'cobird:geometry:update_manifest:badManifest', ...
        'The geometry manifest has an unsupported format.');
else
    manifest = struct('schema', 'COBIRD_GEOMETRY_MANIFEST', ...
        'generated_utc', '', 'models', struct([]));
end

obj_path = fullfile(model_root, [model_name, '.obj']);
[folder, stem] = fileparts(obj_path);
[~, model_id] = fileparts(folder);
[vertices, faces] = cobird.geometry.load_obj_mesh(obj_path);
reabsorb_path = fullfile(folder, [stem, '_reabsorb_info.mat']);
shadow_path = fullfile(folder, [stem, '_shadow_horizon.mat']);
view_factor_path = fullfile(folder, [stem, '_VF.mat']);
assert(all(cellfun(@isfile, {reabsorb_path, shadow_path, view_factor_path})), ...
    'cobird:geometry:update_manifest:missingProducts', ...
    'All three geometry caches must exist before recording model %s.', model_name);

shadow = load(shadow_path, 'schema', 'ray_offset_m', 'N_az', 'N_el');
assert(isfield(shadow, 'schema') && strcmp(shadow.schema, 'COBIRD_SHADOW_HORIZON'), ...
    'cobird:geometry:update_manifest:badShadowSchema', ...
    'The shadow cache does not use the supported format: %s', shadow_path);
view_data = load(view_factor_path, 'VF');
assert(isfield(view_data, 'VF') && isfield(view_data.VF, 'params') && ...
    isfield(view_data.VF, 'schema') && strcmp(view_data.VF.schema, 'M2_VIEW_FACTOR_EXCHANGE'), ...
    'cobird:geometry:update_manifest:badViewFactorSchema', ...
    'The view-factor cache does not use the supported format: %s', view_factor_path);

item = struct();
item.model_id = model_id;
item.modelname = local_relative_file(model_root, fullfile(folder, stem));
item.obj_file = local_relative_file(model_root, obj_path);
item.obj_sha256 = cobird.io.compute_file_sha256(obj_path);
item.vertex_count = size(vertices, 1);
item.face_count = size(faces, 1);
item.reabsorb_file = local_relative_file(model_root, reabsorb_path);
item.reabsorb_sha256 = cobird.io.compute_file_sha256(reabsorb_path);
item.shadow_file = local_relative_file(model_root, shadow_path);
item.shadow_sha256 = cobird.io.compute_file_sha256(shadow_path);
item.shadow_N_az = shadow.N_az;
item.shadow_N_el = shadow.N_el;
item.shadow_schema = shadow.schema;
item.shadow_ray_offset_m = shadow.ray_offset_m;
item.view_factor_file = local_relative_file(model_root, view_factor_path);
item.view_factor_sha256 = cobird.io.compute_file_sha256(view_factor_path);
item.view_factor_method = view_data.VF.method;
item.view_factor_F_threshold = view_data.VF.params.F_threshold;
item.view_factor_ray_offset_m = view_data.VF.params.ray_offset_m;

if isempty(manifest.models)
    manifest.models = item;
else
    matches = strcmp({manifest.models.modelname}, item.modelname);
    assert(nnz(matches) <= 1, 'cobird:geometry:update_manifest:duplicateModel', ...
        'The manifest contains duplicate records for model %s.', model_name);
    item = orderfields(item, manifest.models);
    if any(matches)
        manifest.models(matches) = item;
    else
        manifest.models(end+1, 1) = item;
    end
end
manifest.generated_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));

manifest_dir = fileparts(manifest_path);
if ~isfolder(manifest_dir)
    mkdir(manifest_dir);
end
file_id = fopen(manifest_path, 'w');
assert(file_id >= 0, 'cobird:geometry:update_manifest:writeFailed', ...
    'Unable to write the geometry manifest: %s', manifest_path);
cleanup = onCleanup(@() fclose(file_id));
fwrite(file_id, jsonencode(manifest, PrettyPrint=true), 'char');
end

function relative = local_relative_file(root, path)
prefix = [char(java.io.File(root).getCanonicalPath()), filesep];
canonical = char(java.io.File(path).getCanonicalPath());
assert(startsWith(canonical, prefix, 'IgnoreCase', ispc), ...
    'cobird:geometry:update_manifest:pathOutsideRoot', 'File is outside the model root: %s', path);
relative = strrep(canonical(numel(prefix)+1:end), '\', '/');
end
