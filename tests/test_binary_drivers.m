function tests = test_binary_drivers
%TEST_BINARY_DRIVERS Exercise physical drivers with a small generated binary.
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'examples', 'minimal_binary'));
folder = tempname;
mkdir(folder);
test_case.TestData.folder = folder;
crater_options = struct('N_az_mesh', 8, 'N_radial', 3, ...
    'N_az_shadow', 360, 'N_el_shadow', 180, 'ray_offset', 1e-6);
cobird.geometry.build_crater_model(folder, 'crater', crater_options);
test_case.TestData.smooth = local_problem(folder, "smooth");
test_case.TestData.rough = local_problem(folder, "rough");
end

function teardownOnce(test_case)
% Only the temporary directory created by this fixture is removed.
rmdir(test_case.TestData.folder, 's');
end

function testFourPhysicalModes(test_case)
for mode = ["smooth", "rough"]
    problem = test_case.TestData.(mode);
    thermal = cobird.run_thermal(problem, time_step_s=10, step_count=3);
    coupled = cobird.run_coupled(problem, time_step_s=10, step_count=3);
    verifyEqual(test_case, thermal.time_s, [0; 10; 20; 30]);
    verifyEqual(test_case, coupled.force_time_s, [0; 10; 20]);
    verifyTrue(test_case, all(isfinite(coupled.final_state)));
    verifyFalse(test_case, isfield(thermal, 'final_state'));
    verifyGreaterThan(test_case, thermal.diagnostics.maximum_mutual_flux_W_m2, 0);
    for body = ["primary", "secondary"]
        value = coupled.final_temperature.(body).smooth;
        verifyTrue(test_case, all(isfinite(value), 'all'));
        verifyGreaterThan(test_case, min(value, [], 'all'), 0);
        if mode == "rough"
            verifyNotEmpty(test_case, coupled.final_temperature.(body).rough);
        end
    end
end
end

function testPhysicalSaveStride(test_case)
problem = test_case.TestData.rough;
dense = cobird.run_coupled(problem, time_step_s=10, step_count=4);
sparse = cobird.run_coupled(problem, time_step_s=10, step_count=4, save_stride=3);
verifyEqual(test_case, dense.final_state, sparse.final_state);
verifyEqual(test_case, dense.final_temperature, sparse.final_temperature);
verifyEqual(test_case, sparse.time_s, [0; 30; 40]);
end

function testRestartUsesActualCraterFile(test_case)
problem = test_case.TestData.rough;
file = string(fullfile(test_case.TestData.folder, 'restart.mat'));
first = cobird.run_coupled(problem, time_step_s=10, step_count=2, checkpoint_file=file);
second = cobird.run_coupled(problem, time_step_s=10, step_count=2, resume_file=file);
whole = cobird.run_coupled(problem, time_step_s=10, step_count=4);
verifyEqual(test_case, first.time_s(end), second.time_s(1));
verifyEqual(test_case, second.final_state, whole.final_state);
verifyEqual(test_case, second.final_temperature, whole.final_temperature);
loaded = load(file, 'snapshot');
verifyEqual(test_case, loaded.snapshot.contract.extras.crater_file_sha256, ...
    cobird.io.compute_file_sha256(fullfile(test_case.TestData.folder, 'crater.mat')));
verifyError(test_case, @() cobird.run_coupled(problem, time_step_s=10, ...
    step_count=0, checkpoint_file=file), 'cobird:runtime:integrate_problem:checkpointExists');
alternate = fullfile(test_case.TestData.folder, 'changed');
mkdir(alternate);
loaded_crater = load(fullfile(test_case.TestData.folder, 'crater.mat'), 'crater');
crater = loaded_crater.crater;
marker = 'different file';
save(fullfile(alternate, 'crater.mat'), 'crater', 'marker');
changed = local_problem(alternate, "rough");
verifyError(test_case, @() cobird.run_coupled(changed, time_step_s=10, ...
    step_count=1, resume_file=file), ...
    'cobird:runtime:validate_resume_contract:craterFingerprintMismatch');
end

function testThermalRestart(test_case)
problem = test_case.TestData.smooth;
file = string(fullfile(test_case.TestData.folder, 'thermal_restart.mat'));
cobird.run_thermal(problem, time_step_s=10, step_count=2, checkpoint_file=file);
continued = cobird.run_thermal(problem, time_step_s=10, step_count=2, resume_file=file);
whole = cobird.run_thermal(problem, time_step_s=10, step_count=4);
verifyEqual(test_case, continued.final_temperature, whole.final_temperature);
verifyError(test_case, @() cobird.run_thermal(problem, time_step_s=20, ...
    step_count=1, resume_file=file), 'cobird:runtime:validate_resume_contract:configMismatch');
end

function testThreeDimensionalDriver(test_case)
problem = local_problem(test_case.TestData.folder, "smooth", "3d");
result = cobird.run_coupled(problem, time_step_s=10, step_count=2);
verifySize(test_case, result.final_state, [13, 1]);
verifyTrue(test_case, all(isfinite(result.final_state)));
verifyEqual(test_case, norm(result.final_state(7:10)), 1, 'AbsTol', 1e-12);
end

function testThreeDimensionalEphemerisVectors(test_case)
environment = struct('sun_position_m', [1.495978707e11, 0, 0], ...
    'EMB_position_m', [4e10, 3e10, 1e10], ...
    'Venus_position_m', [-6e10, 2e10, 1e10], ...
    'Jupiter_position_m', [5e11, -3e11, 2e10], 'primary_to_inertial', eye(3));
settings = struct('environment', @(~) environment, ...
    'forces', struct('Sun', true, 'EMB', true, 'Venus', true, 'Jupiter', true));
row_problem = local_problem(test_case.TestData.folder, "smooth", "3d", false, settings);
fields = {'sun_position_m', 'EMB_position_m', 'Venus_position_m', 'Jupiter_position_m'};
for index = 1:numel(fields)
    field = fields{index};
    environment.(field) = environment.(field)';
end
settings.environment = @(~) environment;
column_problem = local_problem(test_case.TestData.folder, "smooth", "3d", false, settings);
radiation = row_problem.evaluate_radiation(0, row_problem.initial_state, row_problem.initial_temperature);
for index = 1:numel(fields)
    verifyEqual(test_case, radiation.environment.(fields{index}), environment.(fields{index}));
end
rows = cobird.run_coupled(row_problem, time_step_s=10, step_count=2);
columns = cobird.run_coupled(column_problem, time_step_s=10, step_count=2);
verifyEqual(test_case, rows.final_state, columns.final_state);
verifyEqual(test_case, rows.final_temperature, columns.final_temperature);
end

function testThreeDimensionalForceOptions(test_case)
for name = ["J2_Dimorphos", "C22_Dimorphos"]
    settings = struct('forces', struct(char(name), false));
    verifyError(test_case, @() local_problem(test_case.TestData.folder, "smooth", "3d", ...
        false, settings), 'cobird:prepare_binary:planarForceOption');
end
verifyError(test_case, @() local_problem(test_case.TestData.folder, "smooth", "planar", ...
    false, struct('forces', struct('secondary_gravity', false))), ...
    'cobird:prepare_binary:threeDimensionalForceOption');
default_problem = local_problem(test_case.TestData.folder, "smooth", "3d");
explicit_problem = local_problem(test_case.TestData.folder, "smooth", "3d", false, ...
    struct('forces', struct('secondary_gravity', true)));
default_result = cobird.run_coupled(default_problem, time_step_s=10, step_count=2);
explicit_result = cobird.run_coupled(explicit_problem, time_step_s=10, step_count=2);
verifyEqual(test_case, default_result.final_state, explicit_result.final_state);
verifyEqual(test_case, default_problem.resume_contract, explicit_problem.resume_contract);
end

function testThreeDimensionalSecondaryGravity(test_case)
settings = struct('forces', struct('Yarkovsky', false, 'secondary_gravity', false));
[problem, params, forces] = local_problem(test_case.TestData.folder, "smooth", "3d", false, settings);
state = problem.initial_state;
state(7:10) = [cos(pi/12); 0; 0; sin(pi/12)];
params.et0 = 0;
inputs = {zeros(3,1), zeros(3,1), zeros(3,1), [0; 0; 1], [], [], [], []};
without = cobird.dynamics.evaluate_3d_coupled_derivative(0, state, params, forces, inputs{:});
forces.secondary_gravity = true;
with = cobird.dynamics.evaluate_3d_coupled_derivative(0, state, params, forces, inputs{:});
default_forces = rmfield(forces, 'secondary_gravity');
default_rhs = cobird.dynamics.evaluate_3d_coupled_derivative(0, state, params, default_forces, inputs{:});
verifyEqual(test_case, default_rhs, with);

r = state(1:3);
w = state(11:13);
rotation = cobird.attitude.convert_quaternion_to_dcm(state(7:10));
verifyEqual(test_case, without(4:6), -params.mu_system*r/norm(r)^3, 'AbsTol', 1e-16);
verifyEqual(test_case, without(11:13), -(params.I \ cross(w, params.I*w)), 'AbsTol', 1e-20);
verifyGreaterThan(test_case, norm(without(11:13)), 0);

% The force must be the negative gradient of the quadrupole potential.
potential = @(position) params.mu_system/(2*norm(position)^3) * ...
    (3*(rotation'*position)'*params.I*(rotation'*position)/norm(position)^2 - trace(params.I));
increment = 1e-4*norm(r);
expected_acceleration = zeros(3,1);
for axis = 1:3
    offset = zeros(3,1);
    offset(axis) = increment;
    expected_acceleration(axis) = -(potential(r+offset)-potential(r-offset))/(2*increment);
end
acceleration = with(4:6)-without(4:6);
verifyEqual(test_case, acceleration, expected_acceleration, 'RelTol', 1e-6, 'AbsTol', 1e-16);

% Gravitational spin and relative-orbit torques must balance.
mass_primary = params.mu_Didymos/params.G_const;
mass_secondary = params.mu_Dimorphos/params.G_const;
reduced_mass = mass_primary*mass_secondary/(mass_primary+mass_secondary);
spin_torque = mass_secondary*rotation*params.I*(with(11:13)-without(11:13));
verifyGreaterThan(test_case, norm(spin_torque), 0);
verifyEqual(test_case, reduced_mass*cross(r, acceleration)+spin_torque, zeros(3,1), ...
    'AbsTol', 1e-10*norm(spin_torque));

% Through the public driver, disabled secondary gravity gives a circular Kepler orbit.
result = cobird.run_coupled(problem, time_step_s=10, step_count=2);
radius = state(1);
motion = sqrt(params.mu_system/radius^3);
angle = motion*result.time_s(end);
expected_orbit = [radius*cos(angle); radius*sin(angle); 0; ...
    -radius*motion*sin(angle); radius*motion*cos(angle); 0];
verifyEqual(test_case, result.final_state(1:6), expected_orbit, 'AbsTol', 1e-9);
end

function testTwentyNodeThermalSolvers(test_case)
temperature = reshape(linspace(190, 250, 60), 3, 20);
inputs = {temperature, [1; 0.5; 0], [10; 2; 0], 1200, ...
    0.9 * 5.670374419e-8, 1e-5, 1e-4, 0.1, 16, 1};
[fixed, fixed_flux] = cobird.thermal.calculate_temperature_fixed20_threads_mex(inputs{:});
[general, general_flux] = cobird.thermal.calculate_temperature_tiled_threads_mex(inputs{:});
verifyEqual(test_case, fixed, general, 'AbsTol', 1e-12);
verifyEqual(test_case, fixed_flux, general_flux, 'AbsTol', 1e-12);
end

function testBodySpecificDepthGrids(test_case)
for mode = ["smooth", "rough"]
    problem = local_problem(test_case.TestData.folder, mode, "planar", false, ...
        struct('body_grids_only', true));
    actual = cobird.run_thermal(problem, time_step_s=10, step_count=2);
    expected = cobird.run_thermal(test_case.TestData.(mode), time_step_s=10, step_count=2);
    verifyEqual(test_case, actual.final_temperature, expected.final_temperature);
end
end

function testEquivalentQuaternionScales(test_case)
settings = struct('use_default_solver', true);
unit = local_problem(test_case.TestData.folder, "smooth", "3d", false, settings);
settings.quaternion_scale = 2;
scaled = local_problem(test_case.TestData.folder, "smooth", "3d", false, settings);
verifyEqual(test_case, scaled.initial_state, unit.initial_state);
verifyEqual(test_case, norm(scaled.initial_state(7:10)), 1, 'AbsTol', 1e-14);
expected = cobird.run_coupled(unit, time_step_s=10, step_count=2);
actual = cobird.run_coupled(scaled, time_step_s=10, step_count=2);
verifyEqual(test_case, actual.final_state, expected.final_state);
verifyEqual(test_case, actual.final_temperature, expected.final_temperature);

trajectory = @(t) local_circular_state_3d(t, unit.initial_state(1), unit.initial_state(13));
state_scale = ones(13, 1);
state_scale(7:10) = 2;
unit_prescribed = local_problem(test_case.TestData.folder, "smooth", "3d", false, ...
    struct('state_at_time', trajectory));
scaled_prescribed = local_problem(test_case.TestData.folder, "smooth", "3d", false, ...
    struct('state_at_time', @(t) state_scale.*trajectory(t)));
expected = cobird.run_thermal(unit_prescribed, time_step_s=10, step_count=2);
actual = cobird.run_thermal(scaled_prescribed, time_step_s=10, step_count=2);
verifyEqual(test_case, actual.final_temperature, expected.final_temperature);
end

function testPlanarStateValidation(test_case)
system = make_synthetic_binary();
contact_radius = system.primary.radius_m+system.secondary.radius_m;
radii = [-system.separation_m, 0, contact_radius/2, contact_radius];
identifiers = ["cobird:prepare_binary:radius", "cobird:prepare_binary:radius", ...
    "cobird:prepare_binary:separation", "cobird:prepare_binary:separation"];
for index = 1:numel(radii)
    state = test_case.TestData.smooth.initial_state;
    state(1) = radii(index);
    verifyError(test_case, @() local_problem(test_case.TestData.folder, "smooth", "planar", ...
        false, struct('initial_state', state)), char(identifiers(index)));
    problem = local_problem(test_case.TestData.folder, "smooth", "planar", false, ...
        struct('state_at_time', @(~) state));
    verifyError(test_case, @() cobird.run_thermal(problem, time_step_s=10, step_count=1), ...
        char(identifiers(index)));
end
end

function testThreeDimensionalQuaternionValidation(test_case)
verifyError(test_case, @() local_problem(test_case.TestData.folder, "smooth", "3d", ...
    false, struct('quaternion_scale', 0)), 'cobird:prepare_binary:quaternion');
valid = local_problem(test_case.TestData.folder, "smooth", "3d");
state = valid.initial_state;
state(7:10) = 0;
problem = local_problem(test_case.TestData.folder, "smooth", "3d", false, ...
    struct('state_at_time', @(~) state));
verifyError(test_case, @() cobird.run_thermal(problem, time_step_s=10, step_count=1), ...
    'cobird:prepare_binary:quaternion');
end

function testDefaultOrbitSolverRestart(test_case)
problem = local_problem(test_case.TestData.folder, "smooth", "3d", false, ...
    struct('use_default_solver', true));
file = string(fullfile(test_case.TestData.folder, 'default_solver_restart.mat'));
first = cobird.run_coupled(problem, time_step_s=10, step_count=1, ...
    start_time_s=100, checkpoint_file=file);
continued = cobird.run_coupled(problem, time_step_s=10, step_count=1, resume_file=file);
whole = cobird.run_coupled(problem, time_step_s=10, step_count=2, start_time_s=100);
verifyEqual(test_case, first.time_s, [100; 110]);
verifyEqual(test_case, continued.time_s, [110; 120]);
verifyEqual(test_case, continued.final_state, whole.final_state);
verifyEqual(test_case, continued.final_temperature, whole.final_temperature);
verifyTrue(test_case, all(isfinite(whole.final_state)));
end

function testThreeDimensionalPrescribedThermal(test_case)
problem = local_problem(test_case.TestData.folder, "smooth", "3d");
result = cobird.run_thermal(problem, time_step_s=10, step_count=2, start_time_s=100);
verifyEqual(test_case, result.time_s, [100; 110; 120]);
verifyEqual(test_case, result.force_time_s, [100; 110]);
verifyFalse(test_case, isfield(result, 'final_state'));
verifyTrue(test_case, all(isfinite(result.final_temperature.secondary.smooth), 'all'));
verifyGreaterThan(test_case, norm(result.final_temperature.secondary.smooth - ...
    problem.initial_temperature.secondary.smooth, 'fro'), 0);
end

function testSeparateTidalParameters(test_case)
names = {'RUN_TIDAL_K2_OVER_Q_PRIMARY', 'RUN_TIDAL_K2_OVER_Q_SECONDARY'};
previous = cellfun(@getenv, names, 'UniformOutput', false);
cleanup = onCleanup(@() local_restore_environment(names, previous));
setenv(names{1}, '');
setenv(names{2}, '');
config = struct('tidal', struct('k2_over_Q', 1e-5, 'k2_over_Q_primary', 2e-5));
[~, ~, primary, secondary] = cobird.runtime.resolve_separate_tidal_parameters(config, struct());
verifyEqual(test_case, [primary, secondary], [2e-5, 1e-5]);
setenv(names{2}, '3e-5');
[~, provenance, primary, secondary] = cobird.runtime.resolve_separate_tidal_parameters(config, struct());
verifyEqual(test_case, [primary, secondary], [2e-5, 3e-5]);
verifyEqual(test_case, provenance.effective_config.tidal.k2_over_Q_secondary, 3e-5);
setenv(names{2}, '0');
verifyError(test_case, @() cobird.runtime.resolve_separate_tidal_parameters(config, struct()), ...
    'cobird:runtime:resolve_separate_tidal_parameters:invalidValue');
end

function testIndependentTidalDynamics(test_case)
names = {'RUN_TIDAL_K2_OVER_Q', 'RUN_TIDAL_K2_OVER_Q_PRIMARY', ...
    'RUN_TIDAL_K2_OVER_Q_SECONDARY', 'RUN_TIDAL_COEFF_A', 'RUN_TIDAL_COEFF_B', ...
    'RUN_TIDAL_I_BZ', 'RUN_TIDAL_RHO_SECONDARY', 'RUN_TIDAL_M_RED'};
previous = cellfun(@getenv, names, 'UniformOutput', false);
cleanup = onCleanup(@() local_restore_environment(names, previous));
for index = 1:numel(names)
    setenv(names{index}, '');
end
values = [1e-5, 1e-5; 2e-5, 1e-5; 1e-5, 2e-5];
results = cell(3, 1);
for index = 1:3
    problem = local_problem(test_case.TestData.folder, "smooth", "planar", false, ...
        struct('tidal_values', values(index, :)));
    results{index} = cobird.run_coupled(problem, time_step_s=1, step_count=1);
end
primary_change = results{2}.final_state-results{1}.final_state;
secondary_change = results{3}.final_state-results{1}.final_state;
% Primary tides transfer spin to the orbit without directly spinning up the secondary.
verifyGreaterThan(test_case, primary_change(4), 0);
verifyEqual(test_case, primary_change(4)+primary_change(6), 0, 'AbsTol', 1e-18);
% A faster-than-synchronous secondary loses spin while its orbit gains angular speed.
verifyGreaterThan(test_case, secondary_change(4), 0);
verifyLessThan(test_case, secondary_change(4)+secondary_change(6), 0);
end

function testTerrainAndNonuniformColumns(test_case)
problem = local_problem(test_case.TestData.folder, "rough", "planar", true);
source = problem.evaluate_radiation(0, problem.initial_state, problem.initial_temperature);
verifyGreaterThan(test_case, min(source.primary.extra), 0);
verifyGreaterThan(test_case, max(source.primary.extra_rough, [], 'all'), 0);
verifyEqual(test_case, source.params.Gamma, 300*2^0.5, 'AbsTol', 1e-12);
result = cobird.run_coupled(problem, time_step_s=10, step_count=2);
verifyTrue(test_case, all(isfinite(result.final_temperature.secondary.rough), 'all'));
end

function local_restore_environment(names, values)
for index = 1:numel(names)
    setenv(names{index}, values{index});
end
end

function [problem, params, forces] = local_problem(folder, surface_mode, geometry, terrain, settings)
if nargin < 3
    geometry = "planar";
end
if nargin < 4
    terrain = false;
end
if nargin < 5
    settings = struct();
end
system = make_synthetic_binary();
primary = system.primary;
secondary = system.secondary;
params = struct('A', 0.1, 'epsilon', 0.9, 'F_s', 1361, 'sigma', 5.670374419e-8, ...
    'rho', 1800, 'Cp', 680, 'Gamma', 300, 'dz', 0.05, 'Nz', 6, 'c', 299792458, ...
    'roughness_f_R', 0.5, 'G_const', 6.67430e-11, ...
    'a2', secondary.axes_m(1)^2, 'b2', secondary.axes_m(2)^2, ...
    'J2A_R2A', 0, 'J2R2_Dim', 0, 'C22R2', 0);
params.k = params.Gamma^2/(params.rho*params.Cp);
params.mu_Didymos = params.G_const*primary.mass_kg;
params.mu_Dimorphos = params.G_const*secondary.mass_kg;
params.mu_system = params.mu_Didymos+params.mu_Dimorphos;
motion = sqrt(params.mu_system/system.separation_m^3);
initial = [system.separation_m; 0; 0; motion; 0; 0];
trajectory = @(t) [system.separation_m; 0; motion*t; motion; 0; 0];
forces = struct('J2_Didymos', false, 'J2_Dimorphos', false, 'C22_Dimorphos', false, ...
    'DART_drag', false, 'tidal', false, 'Yarkovsky', true, ...
    'Sun', false, 'EMB', false, 'Venus', false, 'Jupiter', false);
if geometry == "3d"
    forces = rmfield(forces, {'J2_Dimorphos', 'C22_Dimorphos'});
    params.I = diag([secondary.axes_m(2)^2+secondary.axes_m(3)^2, ...
        secondary.axes_m(1)^2+secondary.axes_m(3)^2, params.a2+params.b2])/5;
    params.GM_Sun = 1.32712440018e20;
    initial = [system.separation_m; 0; 0; 0; motion*system.separation_m; 0; 1; 0; 0; 0; 1e-6; 2e-6; motion];
    trajectory = @(t) local_circular_state_3d(t, system.separation_m, motion);
    if isfield(settings, 'quaternion_scale')
        initial(7:10) = settings.quaternion_scale*initial(7:10);
    end
end
params.Gamma_0 = 300;
params.Gamma_r0_AU = 2;
params.Gamma_xi = -0.5;
if terrain
    primary = local_cavity(primary);
    secondary = local_cavity(secondary);
    secondary.dz = [0, 0.02, 0.05, 0.1, 0.18, 0.3];
    secondary.Nz = 6;
end
if isfield(settings, 'body_grids_only') && settings.body_grids_only
    primary.dz = params.dz;
    primary.Nz = params.Nz;
    secondary.dz = params.dz;
    secondary.Nz = params.Nz;
    params = rmfield(params, {'dz', 'Nz'});
end
if isfield(settings, 'tidal_values')
    forces.tidal = true;
    forces.Yarkovsky = false;
    initial(6) = 1e-4;
    params.tidal = struct('k2_over_Q', 1e-5, ...
        'k2_over_Q_primary', settings.tidal_values(1), ...
        'k2_over_Q_secondary', settings.tidal_values(2), ...
        'coeff_A', 1.5*params.G_const*secondary.mass_kg^2*primary.radius_m^5, ...
        'coeff_B', 1.5*params.G_const*primary.mass_kg^2*secondary.radius_m^5, ...
        'I_Bz', secondary.mass_kg*(params.a2+params.b2)/5, ...
        'rho_secondary', secondary.density_kg_m3, ...
        'm_red', primary.mass_kg*secondary.mass_kg/(primary.mass_kg+secondary.mass_kg));
end
if isfield(settings, 'initial_state')
    initial = settings.initial_state;
end
if isfield(settings, 'state_at_time')
    trajectory = settings.state_at_time;
end
if isfield(settings, 'forces')
    names = fieldnames(settings.forces);
    for index = 1:numel(names)
        forces.(names{index}) = settings.forces.(names{index});
    end
end
options = struct('initial_state', initial, 'initial_temperature_K', 220, 'geometry', geometry, ...
    'environment', @(~) struct('sun_position_m', [1.495978707e11; 0; 0], ...
        'primary_to_inertial', eye(3)), 'environment_id', "fixed_1au_sun_and_attitude", ...
    'state_at_time', trajectory, 'orbit_id', "circular_test_orbit", 'forces', forces, ...
    'surface_mode', surface_mode, 'terrain_heating', terrain, 'reverse_heating', true, ...
    'variable_thermal_inertia', terrain, 'crater_file', string(fullfile(folder, 'crater.mat')));
if isfield(settings, 'environment')
    options.environment = settings.environment;
end
if ~isfield(settings, 'use_default_solver') || ~settings.use_default_solver
    options.orbit_solver = "rk4";
end
arguments_cell = namedargs2cell(options);
problem = cobird.prepare_binary(primary, secondary, params, arguments_cell{:});
end

function state = local_circular_state_3d(time_s, radius, motion)
angle = motion*time_s;
state = [radius*cos(angle); radius*sin(angle); 0; ...
    -radius*motion*sin(angle); radius*motion*cos(angle); 0; ...
    cos(angle/2); 0; 0; sin(angle/2); 0; 0; motion];
end

function body = local_cavity(body)
% Two mutually facing triangles provide a nonzero M2 terrain path.
body.vertices = [-1, 0, 0; -1, 1, 0; -1, 0, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0];
body.faces = [1, 2, 3; 4, 5, 6];
[body.normal, body.area, body.centroid] = ...
    cobird.geometry.compute_face_geometry(body.vertices, body.faces);
body.radius_m = 2;
body.shadow.occ_rev_map = zeros(2, 1, 'uint32');
body.reabsorb = struct('has_candidates', false(2, 1));
body.occultation = struct('no_occ', true(2, 1), 'occ_idx', zeros(0, 1), ...
    'x_local', zeros(0, 3), 'y_local', zeros(0, 3), 'Iv_mat', zeros(0, 3));
body.view_factor = struct('schema', 'M2_VIEW_FACTOR_EXCHANGE', 'N_faces', 2, ...
    'area', body.area, 'pair_i', 1, 'pair_j', 2, 'exchange_area', 1/(16*pi));
end
