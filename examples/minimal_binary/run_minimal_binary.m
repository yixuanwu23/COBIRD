function result = run_minimal_binary(options)
%RUN_MINIMAL_BINARY Run a deterministic synthetic coupled-binary example.
% The example uses generated meshes, an analytic Sun direction, explicit
% heat integration, directional rough radiation, photon recoil, and planar
% orbit-libration dynamics. It requires no external data or MICE.

arguments
    options.step_count (1,1) double {mustBeInteger, mustBePositive} = 48
    options.time_step_s (1,1) double {mustBeFinite, mustBePositive} = 300
    options.verbose (1,1) logical = true
    options.coupled (1,1) logical = true
    options.surface_mode (1,1) string {mustBeMember(options.surface_mode,["smooth","rough"])} = "rough"
end

example_root = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(example_root));
addpath(fullfile(project_root, 'src'));

system = make_synthetic_binary();
primary = system.primary;
secondary = system.secondary;
params = local_parameters(primary, secondary);
forces = local_force_switches();

cache_root = tempname;
mkdir(cache_root);
cache_cleanup = onCleanup(@() local_remove_directory(cache_root));
crater_options = struct( ...
    'N_az_mesh', 8, ...
    'N_radial', 3, ...
    'N_az_shadow', 360, ...
    'N_el_shadow', 180, ...
    'ray_offset', 1e-6, ...
    'opening_area_mode', 'mesh_aperture');
if options.surface_mode == "rough"
    cobird.geometry.build_crater_model( ...
        cache_root, 'minimal_roughness_crater', crater_options);
end
separation = system.separation_m;
mean_motion = sqrt(params.mu_system/separation^3);
initial_state = [separation;0;0;mean_motion;0;0];
problem = cobird.prepare_binary(primary,secondary,params, ...
    initial_state=initial_state, initial_temperature_K=220, ...
    environment=@(t) local_environment(t,2*pi/(6*3600)), ...
    environment_id="synthetic_1au_6h_spin", ...
    state_at_time=@(t) [separation;0;mean_motion*t;mean_motion;0;0], ...
    orbit_id="synthetic_circular_orbit", forces=forces, ...
    terrain_heating=false, reverse_heating=true, ...
    surface_mode=options.surface_mode, ...
    crater_file=string(fullfile(cache_root,'minimal_roughness_crater.mat')), ...
    thermal_solver="explicit", orbit_solver="rk4", ...
    orbit_substep_s=options.time_step_s);
if options.coupled
    solution = cobird.run_coupled(problem, time_step_s=options.time_step_s, step_count=options.step_count);
    state = solution.final_state;
    state_history = cell2mat(cellfun(@transpose,solution.state,'UniformOutput',false));
else
    solution = cobird.run_thermal(problem, time_step_s=options.time_step_s, step_count=options.step_count);
    state_history = [repmat(separation,numel(solution.time_s),1),zeros(numel(solution.time_s),1), ...
        mean_motion*solution.time_s,repmat(mean_motion,numel(solution.time_s),1),zeros(numel(solution.time_s),2)];
    state = state_history(end,:)';
end
temperature_primary = solution.final_temperature.primary.smooth;
temperature_secondary = solution.final_temperature.secondary.smooth;
all_temperatures = [temperature_primary(:);temperature_secondary(:); ...
    solution.final_temperature.primary.rough(:);solution.final_temperature.secondary.rough(:)];
maximum_shadowed_faces = solution.diagnostics.maximum_shadowed_faces;
maximum_mutual_flux = solution.diagnostics.maximum_mutual_flux_W_m2;
maximum_visible_faces = solution.diagnostics.maximum_visible_faces;
maximum_recoil = solution.diagnostics.maximum_recoil_N;
maximum_rough_recoil_difference = solution.diagnostics.maximum_rough_recoil_difference_N;
surface_temperatures = [temperature_primary(:,1); temperature_secondary(:,1)];
surface_temperature_range = max(surface_temperatures) - min(surface_temperatures);
checks = struct();
checks.geometry = all(primary.area > 0) && all(secondary.area > 0);
checks.eclipse = maximum_shadowed_faces > 0;
checks.mutual_heating = maximum_visible_faces > 0 && maximum_mutual_flux > 0;
checks.thermal = all(isfinite(all_temperatures)) && min(all_temperatures) > 0 && ...
    max(all_temperatures) < 1000 && surface_temperature_range > 1;
checks.roughness = options.surface_mode == "smooth" || maximum_rough_recoil_difference > 0;
checks.recoil = isfinite(maximum_recoil) && maximum_recoil > 0;
checks.dynamics = all(isfinite(state)) && state(1) > 0 && norm(state - initial_state) > 0;

result = struct();
result.schema = 'COBIRD_MINIMAL_BINARY';
result.mutual_geometry = 'common_centerline';
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
result.maximum_rough_recoil_difference_N = maximum_rough_recoil_difference;
result.checks = checks;
result.solution = solution;

check_values = struct2cell(checks);
if ~all([check_values{:}])
    failed_names = fieldnames(checks);
    failed_names = failed_names(~[check_values{:}]);
    error('cobird:example:minimalBinaryFailed', ...
        'Minimal binary checks failed: %s', strjoin(failed_names, ', '));
end

if options.verbose
    local_print_check('Geometry', checks.geometry);
    local_print_check('Eclipse', checks.eclipse);
    local_print_check('Mutual heating', checks.mutual_heating);
    local_print_check('Thermal evolution', checks.thermal);
    if options.surface_mode == "rough"
        local_print_check('Thermal roughness', checks.roughness);
    end
    local_print_check('Radiative recoil', checks.recoil);
    if options.coupled
        local_print_check('State progression', checks.dynamics);
    else
        local_print_check('Orbit sampling', checks.dynamics);
    end
    fprintf('Temperature range: %.2f--%.2f K\n', ...
        result.minimum_temperature_K, result.maximum_temperature_K);
end
end

function params = local_parameters(primary, secondary)
params = struct();
params.A = 0.1;
params.epsilon = 0.9;
params.F_s = 1361;
params.sigma = 5.670374419e-8;
params.rho = 1800;
params.Cp = 680;
params.Gamma = 300;
params.k = params.Gamma^2 / (params.rho * params.Cp);
params.dz = 0.05;
params.Nz = 6;
params.c = 299792458;
params.roughness_f_R = 0.5;
params.G_const = 6.67430e-11;
params.mu_Didymos = params.G_const * primary.mass_kg;
params.mu_Dimorphos = params.G_const * secondary.mass_kg;
params.mu_system = params.mu_Didymos + params.mu_Dimorphos;
params.a2 = secondary.axes_m(1)^2;
params.b2 = secondary.axes_m(2)^2;
params.J2A_R2A = 0;
params.J2R2_Dim = 0;
params.C22R2 = 0;
end

function forces = local_force_switches()
forces = struct();
forces.J2_Didymos = false;
forces.J2_Dimorphos = false;
forces.C22_Dimorphos = false;
forces.DART_drag = false;
forces.tidal = false;
forces.Yarkovsky = true;
end

function local_print_check(label, passed)
if passed
    status = 'PASS';
else
    status = 'FAIL';
end
fprintf('%-20s %s\n', label, status);
end

function local_remove_directory(path_to_remove)
if isfolder(path_to_remove)
    rmdir(path_to_remove, 's');
end
end

function environment = local_environment(time_s, spin_rate)
angle = spin_rate*time_s;
environment = struct('sun_position_m',1.495978707e11*[1;0;0], ...
    'primary_to_inertial',[cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1]);
end
