function tests = test_geometry_manifest
%TEST_GEOMETRY_MANIFEST Verify independent registration of model caches.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
end

function setup(test_case)
root = tempname;
mkdir(root);
test_case.TestData.root = root;
local_create_model(root, 'primary');
local_create_model(root, 'secondary');
end

function teardown(test_case)
if isfolder(test_case.TestData.root)
    rmdir(test_case.TestData.root, 's');
end
end

function testUpdatePreservesChangedSiblingRecord(test_case)
root = test_case.TestData.root;
first = 'processed/primary/shape';
second = 'processed/secondary/shape';
cobird.geometry.update_manifest(root, first);
manifest = cobird.geometry.update_manifest(root, second);
record = manifest.models(strcmp({manifest.models.modelname}, second));
local_write_mesh(fullfile(root, [second, '.obj']), true);

manifest = cobird.geometry.update_manifest(root, first);
verifyEqual(test_case, manifest.models(strcmp({manifest.models.modelname}, second)), record);
verifyError(test_case, @() cobird.geometry.prepare_model_cache(root, second, 'check'), ...
    'cobird:geometry:prepare_model_cache:objChanged');
cobird.geometry.prepare_model_cache(root, first, 'check', local_settings());
end

function testUpdatePreservesSiblingCacheHashes(test_case)
root = test_case.TestData.root;
first = 'processed/primary/shape';
second = 'processed/secondary/shape';
cobird.geometry.update_manifest(root, first);
cobird.geometry.update_manifest(root, second);
reabsorb_info = struct('face_occlusion_info', {{[0, 0, 1]}}, 'candidate_faces', {{[]}});
save(fullfile(root, [second, '_reabsorb_info.mat']), 'reabsorb_info');
cobird.geometry.update_manifest(root, first);
verifyError(test_case, @() cobird.geometry.prepare_model_cache(root, second, 'check'), ...
    'cobird:geometry:prepare_model_cache:cacheChanged');
end

function testUpdateRegistersOnlySelectedModel(test_case)
root = test_case.TestData.root;
manifest = cobird.geometry.update_manifest(root, 'processed/primary/shape');
verifyEqual(test_case, {manifest.models.modelname}, {'processed/primary/shape'});
verifyError(test_case, @() cobird.geometry.prepare_model_cache( ...
    root, 'processed/secondary/shape', 'check'), ...
    'cobird:geometry:prepare_model_cache:missingManifestEntry');
end

function testRejectsStaleReabsorptionMethod(test_case)
root = test_case.TestData.root;
model = 'processed/primary/shape';
path = fullfile(root, [model, '_reabsorb_info.mat']);
data = load(path, 'reabsorb_info');
reabsorb_info = data.reabsorb_info;
reabsorb_info.method = 'unsupported';
save(path, 'reabsorb_info');
% Register the matching hash so rejection specifically tests the method.
cobird.geometry.update_manifest(root, model);
verifyError(test_case, @() cobird.geometry.prepare_model_cache(root, model, 'check'), ...
    'cobird:geometry:validate_manifest:badReabsorptionMethod');
end

function local_create_model(root, name)
folder = fullfile(root, 'processed', name);
mkdir(folder);
local_write_mesh(fullfile(folder, 'shape.obj'), false);
reabsorb_info = struct('method', 'triangle_horizon', ...
    'face_occlusion_info', {{[]}}, 'candidate_faces', {{[]}});
save(fullfile(folder, 'shape_reabsorb_info.mat'), 'reabsorb_info');
schema = 'COBIRD_SHADOW_HORIZON';
ray_offset_m = 1e-6;
N_az = 8;
N_el = 4;
shadow_table = false(0, N_az*N_el);
occ_rev_map = uint32(0);
save(fullfile(folder, 'shape_shadow_horizon.mat'), ...
    'schema', 'ray_offset_m', 'N_az', 'N_el', 'shadow_table', 'occ_rev_map');
VF = struct('schema', 'M2_VIEW_FACTOR_EXCHANGE', 'method', 'M2_REZAC_ZHAO_2020', ...
    'N_faces', 1, 'area', 0.5, 'exchange_area', [], ...
    'params', struct('F_threshold', 1e-4, 'ray_offset_m', 1e-6));
save(fullfile(folder, 'shape_VF.mat'), 'VF');
end

function local_write_mesh(path, rotate)
file_id = fopen(path, 'w');
assert(file_id >= 0, 'Unable to create the temporary mesh.');
cleanup = onCleanup(@() fclose(file_id));
if rotate
    fprintf(file_id, 'v 0 0 0\nv 0 1 0\nv -1 0 0\nf 1 2 3\n');
else
    fprintf(file_id, 'v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n');
end
end

function settings = local_settings()
settings = struct('N_az', 8, 'N_el', 4, 'F_threshold', 1e-4, 'ray_offset_m', 1e-6);
end
