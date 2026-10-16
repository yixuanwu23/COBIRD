function problem = prepare_binary(primary, secondary, params, options)
%PREPARE_BINARY Assemble the physical model used by both public drivers.
% Use SI units and fixed meshes centred on each body's centre of mass.
% Triangle winding must give outward normals. Each body contains:
%   vertices [Nv,3] [m], faces [Nf,3] triangle indices, normal [Nf,3] unit vectors
%   area [Nf,1] [m^2], centroid [Nf,3] [m]
%   radius_m: enclosing-sphere radius used for separation and eclipse checks
%   shadow: loaded horizon table; reabsorb: face visibility candidates
%   occultation: output of build_occultation_precomputation
%   view_factor: VF from build_view_factor_cache when terrain_heating is true
%   dz, Nz: optional body-specific depth spacing/nodes and column length
% Missing dz or Nz is taken from params. dz is a positive spacing [m], or
% Nz increasing depth nodes [m] starting at zero for the implicit solver.
% All caches must use the same mesh and metre-based geometry; see
% data/shape_models/dart_pre_impact/README.md for their construction.
%
% Required options are initial_state, environment, environment_id and forces.
% The params structure supplies physical constants and material properties:
%   A, epsilon: albedo and emissivity; F_s: solar flux at 1 au [W/m^2]
%   sigma: Stefan-Boltzmann constant [W/(m^2 K^4)]
%   rho: thermal material density [kg/m^3]; Cp: heat capacity [J/(kg K)]
%   k: thermal conductivity [W/(m K)]; c: speed of light [m/s]
%   dz, Nz: depth spacing/nodes and count, unless supplied on each body
%   mu_Didymos, mu_Dimorphos: primary and secondary GM [m^3/s^2]
%   mu_system: their sum; G_const: gravitational constant [m^3/(kg s^2)]
% Material properties apply to both bodies. Field names containing Didymos
% refer to the primary; those containing Dimorphos refer to the secondary.
% With variable_thermal_inertia=true, also supply Gamma_0 [J/(m^2 K s^0.5)],
% Gamma_r0_AU and Gamma_xi. Gamma=Gamma_0*(distance_au/Gamma_r0_AU)^Gamma_xi
% and k=Gamma^2/(rho*Cp) are then evaluated at each thermal step.
%
% Planar state: [r; dr/dt; theta; dtheta/dt; phi; dphi/dt], where theta is
% the primary-to-secondary longitude in orbit_basis. At phi=0, the secondary
% +x axis points toward the primary and +z follows the orbit normal; positive
% phi rotates +x toward +y. orbit_basis columns are the radial direction at
% theta=0, tangential direction and orbit normal. Angles are in radians.
% Planar params additionally require a2, b2 (squared secondary x/y semiaxes),
% J2A_R2A, J2R2_Dim and C22R2 (gravity coefficients times reference radius
% squared). All five fields have units m^2, including when their terms are off.
% Planar forces contains logical J2_Didymos, J2_Dimorphos, C22_Dimorphos,
% DART_drag, tidal and Yarkovsky fields. Yarkovsky enables both recoil force
% and recoil torque in the dynamics.
%
% 3D state: [r(3); v(3); q(4); omega(3)], where r/v are primary-to-secondary
% inertial vectors [m, m/s], q is scalar-first and body-to-inertial, and omega is
% secondary body-frame angular velocity [rad/s]. Nonzero input q is normalized.
% 3D params additionally require I, the secondary body-frame inertia tensor
% divided by its mass [m^2]. Supply dimensionless J2_Didymos and R_Didymos [m]
% when primary J2 is enabled, and GM_Sun [m^3/s^2] when solar gravity is enabled.
% 3D forces contains logical J2_Didymos, DART_drag, tidal, Yarkovsky, Sun, EMB,
% Venus and Jupiter fields. Optional secondary_gravity (default true) enables
% the secondary quadrupole orbital acceleration and gravity-gradient torque
% together. I still supplies the rigid-body inertia when this switch is false.
% The J2_Dimorphos/C22_Dimorphos switches are planar-only; secondary_gravity
% is 3D-only. Yarkovsky enables both radiative force and torque in either mode.
%
% If DART_drag is enabled, supply t_imp, drag_tau and drag_cutoff [s], and
% drag_A [rad/s^2]. t_imp uses the same epoch as the simulation time.
% If tidal is enabled, params.tidal supplies k2_over_Q or both
% k2_over_Q_primary and k2_over_Q_secondary; body values override the shared one.
% Planar tides assume prograde primary rotation faster than the orbital motion.
% Planar tides also require tidal.coeff_A and coeff_B, where
% coeff_A=(3/2)*G*M_secondary^2*R_primary^5 and
% coeff_B=(3/2)*G*M_primary^2*R_secondary^5, plus tidal.I_Bz [kg m^2],
% tidal.rho_secondary [kg/m^3] and tidal.m_red (reduced mass) [kg].
% For 3D tidal parameters, see help cobird.dynamics.compute_tidal_response.
%
% terrain_heating (default true) enables exchange within each body.
% mutual_heating (default true) enables primary-to-secondary infrared heating;
% reverse_heating (default false) adds secondary-to-primary infrared heating.
% surface_mode="rough" requires params.roughness_f_R, an area fraction in [0,1],
% and crater_file, a MAT file containing crater. Create it with
% cobird.geometry.build_crater_model(cache_root, name, options), then supply
% fullfile(cache_root, [name, '.mat']). Both bodies use this crater and fraction.
% initial_temperature_K (default 200) initializes every thermal depth node.
% Alternatively, initial_temperature.primary and .secondary each contain
% smooth [Nf,Nz] and rough [Nc*Nf,Nz] temperatures [K]. Column 1 is the surface;
% rough rows group Nc crater facets within each parent facet. In smooth mode,
% rough is empty. Nc is the number of facets in crater.
%
% environment(t) provides sun_position_m (primary-to-Sun, inertial) and
% primary_to_inertial (a body-to-inertial rotation matrix). Enabled 3D perturbers
% need EMB_position_m, Venus_position_m or Jupiter_position_m in that frame.
% The primary attitude is prescribed; only the secondary attitude is evolved.
% Time is seconds from the caller's common epoch. SPICE-derived callbacks
% must convert that time to ET and positions from kilometres to metres.
% For run_thermal, state_at_time(t) supplies a state in the layout above, and
% orbit_id identifies its data/settings. environment_id similarly identifies
% the environment callback. Update these identifiers when their inputs change.

arguments
    primary (1, 1) struct
    secondary (1, 1) struct
    params (1, 1) struct
    options.initial_state (:, 1) double
    options.environment (1, 1) function_handle
    options.environment_id (1, 1) string
    options.state_at_time = []
    options.orbit_id (1, 1) string = ""
    options.geometry (1, 1) string {mustBeMember(options.geometry, ["planar", "3d"])} = "planar"
    options.orbit_basis (3, 3) double = eye(3)
    options.forces (1, 1) struct
    options.surface_mode (1, 1) string {mustBeMember(options.surface_mode, ["smooth", "rough"])} = "smooth"
    options.terrain_heating (1, 1) logical = true
    options.mutual_heating (1, 1) logical = true
    options.reverse_heating (1, 1) logical = false
    options.variable_thermal_inertia (1, 1) logical = false
    options.crater_file (1, 1) string = ""
    options.initial_temperature (1, 1) struct = struct()
    options.initial_temperature_K (1, 1) double {mustBeFinite, mustBePositive} = 200
    options.thermal_solver (1, 1) string {mustBeMember(options.thermal_solver, ["implicit", "explicit"])} = "implicit"
    options.solver_config (1, 1) struct = struct('thread_count', 1, 'tile_size', 16)
    options.orbit_solver (1, 1) string {mustBeMember(options.orbit_solver, ["ode15s", "rk4"])} = "ode15s"
    options.orbit_substep_s (1, 1) double {mustBeFinite, mustBePositive} = 10
end

assert(strlength(options.environment_id) > 0, ...
    'cobird:prepare_binary:environmentId', ...
    'Provide an environment_id identifying the ephemeris and primary attitude.');
assert(norm(options.orbit_basis' * options.orbit_basis - eye(3), 'fro') < 1e-10 && ...
    det(options.orbit_basis) > 0, 'cobird:prepare_binary:basis', ...
    'orbit_basis must be a right-handed orthonormal frame.');
options.initial_state = local_prepare_state(options.initial_state, options.geometry, primary, secondary);
if options.geometry == "3d"
    assert(~any(isfield(options.forces, {'J2_Dimorphos', 'C22_Dimorphos'})), ...
        'cobird:prepare_binary:planarForceOption', ...
        'For 3D, remove J2_Dimorphos/C22_Dimorphos and use forces.secondary_gravity.');
    if ~isfield(options.forces, 'secondary_gravity')
        options.forces.secondary_gravity = true;
    end
    validateattributes(options.forces.secondary_gravity, {'logical'}, {'scalar'});
else
    assert(~isfield(options.forces, 'secondary_gravity'), ...
        'cobird:prepare_binary:threeDimensionalForceOption', ...
        'For planar dynamics, use J2_Dimorphos/C22_Dimorphos instead of secondary_gravity.');
end
rough = options.surface_mode == "rough";
cobird.thermal.require_mex_solvers( ...
    thermal=options.thermal_solver == "implicit", radiation=rough);
provenance = struct();
if options.forces.tidal
    config = struct('tidal', params.tidal);
    if ~isfield(config.tidal, 'k2_over_Q')
        assert(all(isfield(config.tidal, {'k2_over_Q_primary', 'k2_over_Q_secondary'})), ...
            'cobird:prepare_binary:tides', 'Specify shared k2/Q or both body values.');
        config.tidal.k2_over_Q = config.tidal.k2_over_Q_primary;
    end
    [config, provenance] = cobird.runtime.resolve_run_configuration(config);
    [config, provenance] = cobird.runtime.resolve_separate_tidal_parameters(config, provenance);
    params.tidal = config.tidal;
end
crater = [];
extras = struct('radiation', cobird.thermal.radiation_model_contract());
if rough
    validateattributes(params.roughness_f_R, {'double'}, ...
        {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
    assert(isfile(options.crater_file), 'cobird:prepare_binary:craterFile', ...
        'Rough mode requires a crater MAT file containing crater.');
    extras.crater_file_sha256 = cobird.io.compute_file_sha256(options.crater_file);
    loaded = load(options.crater_file, 'crater');
    crater = loaded.crater;
end
primary = local_prepare_body(primary, params, options, crater);
secondary = local_prepare_body(secondary, params, options, crater);
temperature = options.initial_temperature;
if isempty(fieldnames(temperature))
    temperature.primary = local_initial_temperature(primary, options.initial_temperature_K, crater);
    temperature.secondary = local_initial_temperature(secondary, options.initial_temperature_K, crater);
end
local_validate_temperature(temperature.primary, primary, crater);
local_validate_temperature(temperature.secondary, secondary, crater);

model = struct('primary', primary, 'secondary', secondary, 'params', params, ...
    'options', options, 'crater', crater);
problem = struct('initial_state', options.initial_state, ...
    'initial_temperature', temperature, 'configuration_provenance', provenance);
problem.evaluate_radiation = @(time_s, state, current) ...
    cobird.thermal.evaluate_binary_radiation(model, time_s, state, current);
problem.evaluate_forcing = @evaluate_prescribed;
problem.advance_temperature = @advance_temperature;
problem.advance_state = @(time_s, state, radiation, dt) ...
    cobird.dynamics.advance_binary_state(time_s, state, radiation, dt, model);
problem.summarize_radiation = @(radiation) radiation.diagnostics;
problem.validate_temperature = @validate_temperature;

% Hash the actual model, not file paths or caller-supplied hash strings.
immutable = rmfield(options, {'environment', 'state_at_time', 'initial_state', ...
    'initial_temperature', 'initial_temperature_K', 'crater_file', 'solver_config'});
extras.primary_sha256 = local_body_hash(primary);
extras.secondary_sha256 = local_body_hash(secondary);
problem.resume_contract = cobird.runtime.build_resume_contract( ...
    'binary', cobird.io.compute_value_sha256(params), immutable, extras);

    function radiation = evaluate_prescribed(time_s, current)
        assert(isa(options.state_at_time, 'function_handle') && strlength(options.orbit_id) > 0, ...
            'cobird:prepare_binary:prescribedOrbit', ...
            'run_thermal requires state_at_time and an orbit_id identifying its data.');
        state = local_prepare_state(options.state_at_time(time_s), options.geometry, primary, secondary);
        radiation = problem.evaluate_radiation(time_s, state, current);
    end

    function next = advance_temperature(~, current, radiation, dt)
        next.primary = local_advance_body(current.primary, primary, radiation.primary, ...
            radiation.params, radiation.distance_au, dt, crater, options);
        next.secondary = local_advance_body(current.secondary, secondary, radiation.secondary, ...
            radiation.params, radiation.distance_au, dt, crater, options);
        validate_temperature(next);
    end

    function validate_temperature(current)
        local_validate_temperature(current.primary, primary, crater);
        local_validate_temperature(current.secondary, secondary, crater);
    end
end

function state = local_prepare_state(state, geometry, primary, secondary)
state_count = 6;
if geometry == "3d"
    state_count = 13;
end
validateattributes(state, {'double'}, {'real', 'finite', 'column', 'numel', state_count});
if geometry == "planar"
    assert(state(1) > 0, 'cobird:prepare_binary:radius', 'Orbital radius must be positive.');
    separation = state(1);
else
    separation = norm(state(1:3));
    quaternion = state(7:10);
    scale = max(abs(quaternion));
    assert(scale > 0, 'cobird:prepare_binary:quaternion', 'The quaternion must be nonzero.');
    quaternion = quaternion / scale;
    state(7:10) = quaternion / norm(quaternion);
end
assert(separation > primary.radius_m + secondary.radius_m, ...
    'cobird:prepare_binary:separation', 'Body reference spheres must not overlap.');
end

function body = local_prepare_body(body, params, options, crater)
required = {'vertices', 'faces', 'normal', 'area', 'centroid', 'shadow', 'reabsorb', 'occultation', 'radius_m'};
assert(all(isfield(body, required)), 'cobird:prepare_binary:body', ...
    'Body geometry or visibility caches are incomplete.');
if ~isfield(body, 'dz')
    body.dz = params.dz;
end
if ~isfield(body, 'Nz')
    body.Nz = params.Nz;
end
validateattributes(body.Nz, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(body.dz, {'numeric'}, {'real', 'finite', 'nonnegative', 'vector'});
if isscalar(body.dz)
    assert(body.dz > 0, 'cobird:prepare_binary:depth', 'Depth spacing must be positive.');
else
    assert(numel(body.dz) == body.Nz && body.dz(1) == 0 && all(diff(body.dz) > 0), ...
        'cobird:prepare_binary:depth', 'Depth nodes must increase strictly from zero.');
    assert(options.thermal_solver == "implicit", 'cobird:prepare_binary:nonuniform', ...
        'Nonuniform depth grids require the implicit solver.');
end
if ~isfield(body, 'exchange')
    body.exchange = [];
end
if options.terrain_heating
    assert(isfield(body, 'view_factor'), 'cobird:prepare_binary:exchange', ...
        'Supply body.view_factor from the M2 cache VF structure.');
    view = body.view_factor;
    paths = cobird.geometry.build_rough_macro_paths(body.vertices, body.faces, view);
    exchange_area = sparse([view.pair_i;view.pair_j], [view.pair_j;view.pair_i], ...
        [view.exchange_area;view.exchange_area], view.N_faces, view.N_faces);
    body.exchange = spdiags(1./body.area, 0, view.N_faces, view.N_faces)*exchange_area;
end
body.rough_model = [];
if ~isempty(crater)
    [body.x_local, body.y_local] = cobird.geometry.build_face_tangent_frame( ...
        body.vertices, body.faces, body.normal);
    body.rough_model = cobird.geometry.prepare_rough_radiation( ...
        '', '', body.vertices, body.faces, crater, body.shadow, params.roughness_f_R, params.c, false);
    if options.terrain_heating
        % Reuse the supplied M2 operator without introducing a second cache path.
        body.rough_model.exchange = cobird.geometry.build_rough_exchange_operator( ...
            body.normal, body.x_local, body.y_local, body.area, body.rough_model.table, paths);
        body.rough_model.terrain_heating = true;
    end
end
end

function hash = local_body_hash(body)
names = {'vertices', 'faces', 'normal', 'area', 'centroid', 'shadow', 'reabsorb', ...
    'occultation', 'radius_m', 'dz', 'Nz', 'exchange', 'rough_model'};
payload = struct();
for index = 1:numel(names)
    payload.(names{index}) = body.(names{index});
end
hash = cobird.io.compute_value_sha256(payload);
end

function temperature = local_initial_temperature(body, initial, crater)
temperature = struct('smooth', initial*ones(numel(body.area), body.Nz), 'rough', []);
if ~isempty(crater)
    temperature.rough = initial*ones(numel(body.area)*numel(crater.area), body.Nz);
end
end

function local_validate_temperature(temperature, body, crater)
assert(isequal(size(temperature.smooth), [numel(body.area), body.Nz]), ...
    'cobird:prepare_binary:temperatureSize', 'Smooth temperature dimensions do not match the body.');
if ~isempty(crater)
    assert(isequal(size(temperature.rough), [numel(body.area)*numel(crater.area), body.Nz]), ...
        'cobird:prepare_binary:temperatureSize', 'Rough temperature dimensions do not match the body.');
else
    assert(isempty(temperature.rough), 'cobird:prepare_binary:temperatureSize', ...
        'Smooth mode does not accept rough temperature columns.');
end
values = [temperature.smooth(:);temperature.rough(:)];
assert(all(isfinite(values)) && all(values > 0), ...
    'cobird:prepare_binary:temperature', ...
    'Temperatures must remain positive and finite; check inputs and the time step.');
end

function next = local_advance_body(current, body, heating, params, distance_au, dt, crater, options)
solver_dt = dt;
if options.thermal_solver == "explicit"
    solver_dt = [];
end
next = current;
next.smooth = cobird.thermal.advance_temperature_columns( ...
    [], heating.mu, current.smooth, params, distance_au, solver_dt, heating.extra, ...
    body.dz, body.Nz, options.solver_config);
if ~isempty(crater)
    next.rough = cobird.thermal.advance_rough_temperature( ...
        current.rough, params, distance_au, solver_dt, crater, heating.direct_rough, ...
        heating.extra_rough, heating.emission_rough, body.dz, body.Nz, options.solver_config);
end
if options.thermal_solver == "explicit"
    next.smooth = current.smooth + dt*next.smooth;
    next.rough = current.rough + dt*next.rough;
end
end
