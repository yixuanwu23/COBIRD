function tests = test_core_smoke
%TEST_CORE_SMOKE Fast tests with no external models, kernels, or MICE.
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
test_root = fileparts(mfilename('fullpath'));
project_root = fileparts(test_root);
addpath(fullfile(project_root, 'src'));
test_case.TestData.project_root = project_root;
end

function testQuaternionRoundTrip(test_case)
angle = deg2rad(37);
rotation = [cos(angle), -sin(angle), 0; ...
            sin(angle),  cos(angle), 0; ...
            0,           0,          1];
quaternion = cobird.attitude.convert_dcm_to_quaternion(rotation);
rebuilt = cobird.attitude.convert_quaternion_to_dcm(quaternion);
verifyEqual(test_case, rebuilt, rotation, 'AbsTol', 1e-12);
verifyEqual(test_case, norm(quaternion), 1, 'AbsTol', 1e-14);
end

function testTriangleGeometry(test_case)
vertices = [0, 0, 0; 2, 0, 0; 0, 2, 0];
faces = [1, 2, 3];
[normal, area, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);
verifyEqual(test_case, normal, [0, 0, 1], 'AbsTol', 1e-14);
verifyEqual(test_case, area, 2, 'AbsTol', 1e-14);
verifyEqual(test_case, centroid, [2/3, 2/3, 0], 'AbsTol', 1e-14);
end

function testRayTriangleIntersection(test_case)
vertex_0 = [0, 0, 0];
vertex_1 = [1, 0, 0];
vertex_2 = [0, 1, 0];
hit = cobird.geometry.intersect_triangle_ray( ...
    [0.25, 0.25, 1], [0, 0, -1], vertex_0, vertex_1, vertex_2);
miss = cobird.geometry.intersect_triangle_ray( ...
    [1.25, 1.25, 1], [0, 0, -1], vertex_0, vertex_1, vertex_2);
verifyTrue(test_case, hit);
verifyFalse(test_case, miss);
end

function testObjLoaderWithTextureIndices(test_case)
path = [tempname, '.obj'];
cleanup = onCleanup(@() local_delete_file(path));
file_id = fopen(path, 'w');
assert(file_id >= 0, 'Unable to create the temporary OBJ fixture.');
fprintf(file_id, 'v 0 0 0\n');
fprintf(file_id, 'v 1 0 0\n');
fprintf(file_id, 'v 0 1 0\n');
fprintf(file_id, 'f 1/4/7 2/5/8 3/6/9\n');
fclose(file_id);

[vertices, faces] = cobird.geometry.load_obj_mesh(path);
verifyEqual(test_case, vertices, [0, 0, 0; 1, 0, 0; 0, 1, 0]);
verifyEqual(test_case, faces, [1, 2, 3]);
end

function testObjLoaderWithOptionalVertexFields(test_case)
path = [tempname, '.obj'];
cleanup = onCleanup(@() local_delete_file(path));
file_id = fopen(path, 'w');
assert(file_id >= 0, 'Unable to create the temporary OBJ fixture.');
fprintf(file_id, 'v 0 0 0\nv  1 0 0 1\nv\t0 1 0\n');
fprintf(file_id, 'f\t1 2 3 # triangle\n');
fclose(file_id);

[vertices, faces] = cobird.geometry.load_obj_mesh(path);
verifyEqual(test_case, vertices, [0, 0, 0; 1, 0, 0; 0, 1, 0]);
verifyEqual(test_case, faces, [1, 2, 3]);
end

function testObjLoaderRejectsInvalidVertex(test_case)
path = [tempname, '.obj'];
cleanup = onCleanup(@() local_delete_file(path));
records = {'v 0 0', 'v 0 0 1bad'};
identifiers = {'invalidVertex', 'invalidNumber'};
for index = 1:numel(records)
    file_id = fopen(path, 'w');
    assert(file_id >= 0, 'Unable to create the temporary OBJ fixture.');
    fprintf(file_id, '%s\nv 1 0 0\nv 0 1 0\nf 1 2 3\n', records{index});
    fclose(file_id);
    verifyError(test_case, @() cobird.geometry.load_obj_mesh(path), ...
        ['cobird:geometry:load_obj_mesh:', identifiers{index}]);
end
end

function testTwoBodyAcceleration(test_case)
acceleration = cobird.dynamics.compute_two_body_acceleration([2; 0; 0], 8);
verifyEqual(test_case, acceleration, [-2; 0; 0], 'AbsTol', 1e-14);
end

function testUnoccludedShadowCache(test_case)
cache_root = tempname;
mkdir(cache_root);
cleanup = onCleanup(@() local_remove_directory(cache_root));
vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0];
faces = [1, 2, 3];
[normal, ~, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);
reabsorb_info = struct('face_occlusion_info', {{[]}}, 'candidate_faces', {{[]}});
cache = cobird.geometry.build_shadow_horizon_cache( ...
    vertices, faces, normal, centroid, reabsorb_info, 8, 4, cache_root, 'triangle');
loaded = cobird.geometry.load_shadow_horizon_cache(cache_root, 'triangle', 8, 4);
verifyEqual(test_case, loaded, cache);
verifySize(test_case, cache.shadow_table, [0, 32]);
verifyEqual(test_case, cache.occ_rev_map, uint32(0));
end

function testReflectedFlux(test_case)
flux = cobird.thermal.compute_reflected_light_flux(0.2, [10; 0; 2.5]);
verifyEqual(test_case, flux, [2; 0; 0.5], 'AbsTol', 1e-14);
end

function testCraterApertureAndDirectionGrid(test_case)
cache_root = tempname;
mkdir(cache_root);
cleanup = onCleanup(@() local_remove_directory(cache_root));
options = struct( ...
    'N_az_mesh', 20, ...
    'N_radial', 4, ...
    'N_az_shadow', 20, ...
    'N_el_shadow', 2, ...
    'ray_offset', 1e-6, ...
    'opening_area_mode', 'mesh_aperture');
crater = cobird.geometry.build_crater_model( ...
    cache_root, 'test_crater', options);
verifyEqual(test_case, crater.opening_area, ...
    3.09016994374947, 'AbsTol', 1e-13);
verifyEqual(test_case, crater.schema, 'COBIRD_HEMISPHERICAL_CRATER');
table = cobird.geometry.build_crater_direction_table(crater);
verifyEqual(test_case, ...
    [table.azimuth_count, table.elevation_count], [360, 91]);
verifyLessThan(test_case, table.partition_error, 2e-7);
end

function testCraterDirectionSeam(test_case)
table = struct('schema', 'COBIRD_CRATER_DIRECTION_TABLE', ...
    'azimuth_count', 360, 'elevation_count', 91, ...
    'elevation_step', pi / 180);
directions = [cosd(30), -1e-20, sind(30); ...
    cosd(89.5), -1e-20, sind(89.5); ...
    1e-16, -1e-40, 1; ...
    cosd(30), 0, sind(30); ...
    0, 0, 1];
[indices, coefficients, active] = ...
    cobird.geometry.crater_direction_stencil(table, directions);
verifyTrue(test_case, all(active));
verifyGreaterThanOrEqual(test_case, indices, 1);
verifyLessThanOrEqual(test_case, indices, ...
    table.azimuth_count * table.elevation_count);
verifyGreaterThanOrEqual(test_case, coefficients, 0);
verifyEqual(test_case, sum(coefficients, 1), ones(1, 5), 'AbsTol', 1e-14);
elevations = repelem(0:90, table.azimuth_count);
interpolated = sum(elevations(indices) .* coefficients, 1);
verifyEqual(test_case, interpolated, [30, 89.5, 90, 30, 90], 'AbsTol', 1e-12);
verifyEqual(test_case, indices(:, 1), indices(:, 4));
verifyEqual(test_case, coefficients(:, 1), coefficients(:, 4));
verifyEqual(test_case, indices(:, 3), indices(:, 5));
verifyEqual(test_case, coefficients(:, 3), coefficients(:, 5));
end

function testPlanarRadiativeCoupling(test_case)
params = struct();
params.mu_system = 6;
params.mu_Didymos = 5;
params.mu_Dimorphos = 1;
params.G_const = 1;
params.a2 = 4;
params.b2 = 1;
params.J2A_R2A = 0;
params.J2R2_Dim = 0;
params.C22R2 = 0;
forces = struct('J2_Didymos', false, 'J2_Dimorphos', false, ...
    'C22_Dimorphos', false, 'DART_drag', false, 'tidal', false, ...
    'Yarkovsky', true);
state = [2; 0; 0; sqrt(6/8); 0; 0];
force = [0.2; -0.1];
moment = 0.03;

with_radiation = cobird.dynamics.evaluate_planar_coupled_derivative( ...
    0, state, params, forces, force, moment);
forces.Yarkovsky = false;
without_radiation = cobird.dynamics.evaluate_planar_coupled_derivative( ...
    0, state, params, forces, force, moment);
inertia_z = (params.a2 + params.b2) / 5;
tangential_acceleration = force(2) / state(1);
expected = [0; force(1); 0; tangential_acceleration; 0; ...
    moment / inertia_z - tangential_acceleration];
verifyEqual(test_case, with_radiation - without_radiation, expected, 'AbsTol', 1e-14);
end

function local_delete_file(path)
if isfile(path)
    delete(path);
end
end

function local_remove_directory(path)
if isfolder(path)
    rmdir(path, 's');
end
end
