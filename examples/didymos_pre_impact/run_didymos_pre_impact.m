function result = run_didymos_pre_impact(options)
%RUN_DIDYMOS_PRE_IMPACT Run the official-shape pre-impact quick example.
% The PDS v003/v004 OBJ files are loaded without modification. Vertex
% coordinates are converted from kilometres to metres only in memory.

arguments
    options.step_count (1,1) double {mustBeInteger, mustBePositive} = 24
    options.time_step_s (1,1) double {mustBeFinite, mustBePositive} = 300
    options.verbose (1,1) logical = true
end

example_root = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(example_root));
addpath(fullfile(project_root, 'src'));
data_root = fullfile(project_root, 'data', 'shape_models', 'dart_pre_impact');

[didymos_record, dimorphos_record] = local_source_records();
primary = local_load_pds_body(data_root, didymos_record);
secondary = local_load_pds_body(data_root, dimorphos_record);
params = local_parameters(primary, secondary);
primary.mass_kg = params.mass_primary_kg;
secondary.mass_kg = params.mass_secondary_kg;
forces = local_force_switches();

separation = params.semi_major_axis_m;
mean_motion = 2*pi/params.orbit_period_s;
initial_state = [separation;0;0;mean_motion;0;0];
problem = cobird.prepare_binary(primary,secondary,params, ...
    initial_state=initial_state, initial_temperature_K=params.T_init_K, ...
    environment=@(t) local_environment(t,2*pi/params.primary_rotation_period_s), ...
    environment_id="didymos_quick_1au_2.2600h_spin", ...
    forces=forces, terrain_heating=false, reverse_heating=true, ...
    thermal_solver="explicit", orbit_solver="rk4", ...
    orbit_substep_s=options.time_step_s);
solution = cobird.run_coupled(problem, ...
    time_step_s=options.time_step_s, step_count=options.step_count);
state = solution.final_state;
state_history = cell2mat(cellfun(@transpose,solution.state,'UniformOutput',false));
temperature_primary = solution.final_temperature.primary.smooth;
temperature_secondary = solution.final_temperature.secondary.smooth;
maximum_shadowed_faces = solution.diagnostics.maximum_shadowed_faces;
maximum_mutual_flux = solution.diagnostics.maximum_mutual_flux_W_m2;
maximum_visible_faces = solution.diagnostics.maximum_visible_faces;
maximum_recoil = solution.diagnostics.maximum_recoil_N;
all_temperatures = [temperature_primary(:); temperature_secondary(:)];
surface_temperatures = [temperature_primary(:, 1); temperature_secondary(:, 1)];
surface_temperature_range = max(surface_temperatures) - min(surface_temperatures);

checks = struct();
checks.source_files = primary.source_verified && secondary.source_verified;
checks.geometry = all(primary.area > 0) && all(secondary.area > 0) && ...
    primary.face_count == didymos_record.face_count && ...
    secondary.face_count == dimorphos_record.face_count;
checks.eclipse = maximum_shadowed_faces > 0;
checks.mutual_heating = maximum_visible_faces > 0 && maximum_mutual_flux > 0;
checks.thermal = all(isfinite(all_temperatures)) && min(all_temperatures) > 0 && ...
    max(all_temperatures) < 1000 && surface_temperature_range > 1;
checks.recoil = isfinite(maximum_recoil) && maximum_recoil > 0;
checks.dynamics = all(isfinite(state)) && ...
    state(1) > primary.radius_m + secondary.radius_m && ...
    norm(state - initial_state) > 0;

result = struct();
result.schema = 'COBIRD_DIDYMOS_PRE_IMPACT';
result.mutual_geometry = 'common_centerline';
result.macroscopic_self_shadowing_enabled = false;
result.roughness_enabled = false;
result.step_count = options.step_count;
result.time_step_s = options.time_step_s;
result.elapsed_simulation_s = options.step_count * options.time_step_s;
result.initial_state = initial_state;
result.final_state = state;
result.state_history = state_history;
result.minimum_temperature_K = min(all_temperatures);
result.maximum_temperature_K = max(all_temperatures);
result.surface_temperature_range_K = surface_temperature_range;
result.maximum_shadowed_faces = maximum_shadowed_faces;
result.maximum_visible_faces = maximum_visible_faces;
result.maximum_mutual_flux_W_m2 = maximum_mutual_flux;
result.maximum_recoil_N = maximum_recoil;
result.inferred_system_bulk_density_kg_m3 = params.system_bulk_density_kg_m3;
result.didymos = local_result_metadata(primary);
result.dimorphos = local_result_metadata(secondary);
result.checks = checks;
result.solution = solution;

check_values = struct2cell(checks);
if ~all([check_values{:}])
    failed_names = fieldnames(checks);
    failed_names = failed_names(~[check_values{:}]);
    error('cobird:example:didymosPreImpactFailed', ...
        'Didymos pre-impact checks failed: %s', strjoin(failed_names, ', '));
end

if options.verbose
    fprintf('Official PDS pre-impact shapes.\n');
    fprintf('Macroscopic terrain self-shadowing and terrain heating are disabled.\n');
    local_print_check('Source integrity', checks.source_files);
    local_print_check('Geometry', checks.geometry);
    local_print_check('Eclipse', checks.eclipse);
    local_print_check('Mutual heating', checks.mutual_heating);
    local_print_check('Thermal evolution', checks.thermal);
    local_print_check('Radiative recoil', checks.recoil);
    local_print_check('State progression', checks.dynamics);
    fprintf('Temperature range: %.2f--%.2f K\n', ...
        result.minimum_temperature_K, result.maximum_temperature_K);
end
end

function body = local_load_pds_body(data_root, record)
obj_path = fullfile(data_root, record.obj_file);
xml_path = fullfile(data_root, record.xml_file);
assert(isfile(obj_path) && isfile(xml_path), ...
    'cobird:example:missingPdsShape', ...
    'Required PDS shape product or XML label is missing: %s', record.name);

obj_hash = cobird.io.compute_file_sha256(obj_path);
xml_hash = cobird.io.compute_file_sha256(xml_path);
assert(strcmp(obj_hash, record.obj_sha256) && strcmp(xml_hash, record.xml_sha256), ...
    'cobird:example:pdsShapeHashMismatch', ...
    'PDS source file integrity check failed for %s.', record.name);

[vertices_km, faces] = cobird.geometry.load_obj_mesh(obj_path);
assert(size(vertices_km, 1) == record.vertex_count && ...
    size(faces, 1) == record.face_count, ...
    'cobird:example:pdsShapeCountMismatch', ...
    'PDS vertex or facet count does not match the archived product for %s.', record.name);

% The archived OBJ remains untouched; only this in-memory array uses metres.
vertices = vertices_km * 1000;
[normal, area, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);
[x_local, y_local] = cobird.geometry.build_face_tangent_frame(vertices, faces, normal);
volume_m3 = local_mesh_volume(vertices, faces);
extent_m = max(vertices, [], 1) - min(vertices, [], 1);
face_count = size(faces, 1);

body = struct();
body.name = record.name;
body.product_id = record.product_id;
body.obj_file = record.obj_file;
body.xml_file = record.xml_file;
body.obj_sha256 = obj_hash;
body.xml_sha256 = xml_hash;
body.source_verified = true;
body.vertex_count = size(vertices, 1);
body.face_count = face_count;
body.extent_m = extent_m;
body.axes_m = extent_m / 2;
body.volume_m3 = volume_m3;
body.radius_m = max(vecnorm(vertices, 2, 2));
body.vertices = vertices;
body.faces = faces;
body.normal = normal;
body.area = area;
body.centroid = centroid;
body.x_local = x_local;
body.y_local = y_local;
body.reabsorb = struct('has_candidates', false(face_count, 1));
body.shadow = struct( ...
    'N_az', 1, 'N_el', 1, 'daz', 2*pi, 'del', pi, ...
    'occ_rev_map', zeros(face_count, 1, 'uint32'), ...
    'shadow_table', false(0, 1));
body.occultation = struct( ...
    'no_occ', true(face_count, 1), ...
    'occ_idx', zeros(0, 1), ...
    'x_local', zeros(0, 3), ...
    'y_local', zeros(0, 3), ...
    'Iv_mat', zeros(0, 3));
end

function volume = local_mesh_volume(vertices, faces)
vertex_1 = vertices(faces(:, 1), :);
vertex_2 = vertices(faces(:, 2), :);
vertex_3 = vertices(faces(:, 3), :);
signed_volume = sum(dot(vertex_1, cross(vertex_2, vertex_3, 2), 2)) / 6;
volume = abs(signed_volume);
assert(isfinite(volume) && volume > 0, ...
    'cobird:example:invalidPdsMeshVolume', ...
    'PDS shape mesh has a nonpositive or invalid enclosed volume.');
end

function params = local_parameters(primary, secondary)
params = struct();
params.A = 0.07;
params.epsilon = 0.9;
params.F_s = 1367.5;
params.sigma = 5.670374419e-8;
params.rho = 2400;
params.Cp = 550;
params.Gamma = 320;
params.k = params.Gamma^2 / (params.rho * params.Cp);
params.dz = 0.01;
params.Nz = 6;
params.T_init_K = 200;
params.c = 299792458;
params.d_sun_AU = 1;
params.G_const = 6.67430e-11;
params.semi_major_axis_m = 1189;
params.orbit_period_s = 11.921473 * 3600;
params.primary_rotation_period_s = 2.2600 * 3600;

mean_motion = 2 * pi / params.orbit_period_s;
params.mu_system = mean_motion^2 * params.semi_major_axis_m^3;
total_mass = params.mu_system / params.G_const;
volume_fraction = primary.volume_m3 / (primary.volume_m3 + secondary.volume_m3);
params.mass_primary_kg = total_mass * volume_fraction;
params.mass_secondary_kg = total_mass * (1 - volume_fraction);
params.mu_Didymos = params.G_const * params.mass_primary_kg;
params.mu_Dimorphos = params.G_const * params.mass_secondary_kg;
params.system_bulk_density_kg_m3 = ...
    total_mass / (primary.volume_m3 + secondary.volume_m3);

primary_axes = primary.axes_m;
secondary_axes = secondary.axes_m;
params.a2 = secondary_axes(1)^2;
params.b2 = secondary_axes(2)^2;
params.J2A_R2A = ...
    (primary_axes(1)^2 + primary_axes(2)^2 - 2*primary_axes(3)^2) / 10;
params.J2R2_Dim = ...
    (secondary_axes(1)^2 + secondary_axes(2)^2 - 2*secondary_axes(3)^2) / 10;
params.C22R2 = (secondary_axes(1)^2 - secondary_axes(2)^2) / 20;
end

function forces = local_force_switches()
forces = struct();
forces.J2_Didymos = true;
forces.J2_Dimorphos = true;
forces.C22_Dimorphos = true;
forces.DART_drag = false;
forces.tidal = false;
forces.Yarkovsky = true;
end

function metadata = local_result_metadata(body)
metadata = struct();
metadata.product_id = body.product_id;
metadata.obj_file = body.obj_file;
metadata.xml_file = body.xml_file;
metadata.obj_sha256 = body.obj_sha256;
metadata.xml_sha256 = body.xml_sha256;
metadata.vertex_count = body.vertex_count;
metadata.face_count = body.face_count;
metadata.extent_m = body.extent_m;
metadata.volume_m3 = body.volume_m3;
end

function [didymos, dimorphos] = local_source_records()
didymos = struct( ...
    'name', 'Didymos', ...
    'product_id', 'didymos_g_9309mm_spc_obj_0000n00000_v003', ...
    'obj_file', 'didymos_g_9309mm_spc_obj_0000n00000_v003.obj', ...
    'xml_file', 'didymos_g_9309mm_spc_obj_0000n00000_v003.xml', ...
    'obj_sha256', 'b5349374186d254ca21b7b00ddd9f6dea7e2762ca6bbff13f43151ee448e3f20', ...
    'xml_sha256', 'f09a8a3364c6fc9e9c6df54206b03e550f0da16a896d3d65c6e5d74b0d525564', ...
    'vertex_count', 24578, ...
    'face_count', 49152);
dimorphos = struct( ...
    'name', 'Dimorphos', ...
    'product_id', 'dimorphos_g_1940mm_spc_obj_0000n00000_v004', ...
    'obj_file', 'dimorphos_g_1940mm_spc_obj_0000n00000_v004.obj', ...
    'xml_file', 'dimorphos_g_1940mm_spc_obj_0000n00000_v004.xml', ...
    'obj_sha256', '833eaeac8c49338abe64012a679b1b9695f774f1b5e4842a8e5b9f0af1d28664', ...
    'xml_sha256', '715a4f105d2284ad818b5a3b5eb64087be8f37a6c997fafc31e88724a2c765f6', ...
    'vertex_count', 24578, ...
    'face_count', 49152);
end

function local_print_check(label, passed)
if passed
    status = 'PASS';
else
    status = 'FAIL';
end
fprintf('%-20s %s\n', label, status);
end

function environment = local_environment(time_s, spin_rate)
angle = spin_rate*time_s;
environment = struct('sun_position_m',1.495978707e11*[1;0;0], ...
    'primary_to_inertial',[cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1]);
end
