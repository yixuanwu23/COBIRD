function radiation = evaluate_binary_radiation(model, time_s, state, temperature)
%EVALUATE_BINARY_RADIATION Evaluate heating and recoil at one interval start.
primary = model.primary;
secondary = model.secondary;
options = model.options;
params = model.params;
environment = options.environment(time_s);
position_fields = {'sun_position_m'};
if options.geometry == "3d"
    for name = ["EMB", "Venus", "Jupiter"]
        if options.forces.(name)
            position_fields{end+1} = char(name + "_position_m");
        end
    end
end
for index = 1:numel(position_fields)
    field = position_fields{index};
    position = environment.(field);
    validateattributes(position, {'double'}, {'real', 'finite', 'vector', 'numel', 3}, ...
        'environment', field);
    environment.(field) = position(:);
end
sun = environment.sun_position_m;
rotation_primary = environment.primary_to_inertial;
assert(norm(sun) > 0 && isequal(size(rotation_primary), [3, 3]) && ...
    norm(rotation_primary'*rotation_primary-eye(3), 'fro') < 1e-10 && det(rotation_primary) > 0, ...
    'cobird:thermal:evaluate_binary_radiation:environment', 'The Sun position or primary rotation is invalid.');
if options.geometry == "planar"
    assert(state(1) > 0, 'cobird:thermal:evaluate_binary_radiation:radius', ...
        'Orbital radius must remain positive.');
    basis = options.orbit_basis;
    radial = cos(state(3))*basis(:, 1) + sin(state(3))*basis(:, 2);
    tangential = -sin(state(3))*basis(:, 1) + cos(state(3))*basis(:, 2);
    position = state(1)*radial;
    rotation_secondary = cobird.attitude.compute_planar_body_rotation( ...
        state(3), state(5), basis(:, 1), basis(:, 2), basis(:, 3));
else
    position = state(1:3);
    rotation_secondary = cobird.attitude.convert_quaternion_to_dcm( ...
        state(7:10)/norm(state(7:10)));
end
assert(norm(position) > primary.radius_m+secondary.radius_m, ...
    'cobird:thermal:evaluate_binary_radiation:separation', 'Body reference spheres overlap.');
distance_au = norm(sun)/1.495978707e11;
sun_primary = rotation_primary'*(sun/norm(sun));
sun_secondary = rotation_secondary'*(sun/norm(sun));
position_primary = rotation_primary'*position;
position_secondary = rotation_secondary'*(-position);
if options.variable_thermal_inertia
    params.Gamma = params.Gamma_0*(distance_au/params.Gamma_r0_AU)^params.Gamma_xi;
    params.k = params.Gamma^2/(params.rho*params.Cp);
end
shadow_primary = cobird.geometry.compute_dimorphos_shadow( ...
    primary.centroid, sun_primary, position_primary, secondary.radius_m);
if isfield(params, 'Sun_radius')
    eclipse = cobird.geometry.compute_occultation_factor(-position, sun, ...
        params.Sun_radius, primary.radius_m);
    shadow_secondary = false(numel(secondary.area), 1);
else
    eclipse = 1;
    shadow_secondary = cobird.geometry.compute_dimorphos_shadow( ...
        secondary.centroid, sun_secondary, position_secondary, primary.radius_m);
end
mu_primary = local_incidence(primary, sun_primary, shadow_primary);
mu_secondary = eclipse*local_incidence(secondary, sun_secondary, shadow_secondary);
primary_source = local_sources(primary, temperature.primary, ...
    mu_primary, sun_primary, params, distance_au, model);
secondary_source = local_sources(secondary, temperature.secondary, ...
    mu_secondary, sun_secondary, params, distance_au, model);
rough = options.surface_mode == "rough";
visible_primary = [];
visible_secondary = [];
cosine_primary = [];
cosine_secondary = [];
weights_primary = [];
weights_secondary = [];
endpoint_primary = [];
endpoint_secondary = [];
if options.mutual_heating
    secondary_to_primary = rotation_primary'*rotation_secondary;
    [visible_primary, visible_secondary, cosine_primary, cosine_secondary, distance] = ...
        cobird.geometry.compute_mutual_visibility( ...
        primary.centroid, primary.normal, primary.reabsorb, ...
        (secondary_to_primary*secondary.centroid')'+position_primary', ...
        (secondary_to_primary*secondary.normal')', secondary.reabsorb, ...
        position_primary, primary.radius_m, secondary.radius_m, ...
        primary.shadow, secondary.shadow, position_secondary/norm(position_secondary));
    if ~isempty(visible_primary)
        weights_primary = cosine_primary.*primary.area(visible_primary)/(pi*distance^2);
        weights_secondary = cosine_secondary.*secondary.area(visible_secondary)/(pi*distance^2);
        if rough
            endpoint_primary = cobird.geometry.build_mutual_rough_endpoint( ...
                primary.rough_model, visible_primary, position_primary/norm(position_primary));
            endpoint_secondary = cobird.geometry.build_mutual_rough_endpoint( ...
                secondary.rough_model, visible_secondary, position_secondary/norm(position_secondary));
        end
    end
end
[mutual_primary, mutual_secondary, rough_primary, rough_secondary] = ...
    cobird.thermal.compute_mutual_heating( ...
    primary_source.emission, secondary_source.emission, ...
    primary_source.emission_rough, secondary_source.emission_rough, ...
    visible_primary, visible_secondary, weights_primary, weights_secondary, ...
    cosine_primary, cosine_secondary, endpoint_primary, endpoint_secondary, ...
    params.epsilon, rough, options.reverse_heating);
primary_source.extra = primary_source.extra+mutual_primary;
secondary_source.extra = secondary_source.extra+mutual_secondary;
primary_source.extra_rough = primary_source.extra_rough+rough_primary;
secondary_source.extra_rough = secondary_source.extra_rough+rough_secondary;
[force_primary, moment_primary, ~] = local_recoil(primary, temperature.primary, primary_source, params, rough);
[force_secondary, moment_secondary, rough_difference] = ...
    local_recoil(secondary, temperature.secondary, secondary_source, params, rough);
force_primary_inertial = rotation_primary*force_primary;
radiation = struct('primary', primary_source, 'secondary', secondary_source, 'params', params, ...
    'distance_au', distance_au, 'environment', environment, ...
    'force_primary_body', force_primary, 'moment_primary_body', moment_primary, ...
    'force_primary_inertial', force_primary_inertial, ...
    'force_secondary_body', force_secondary, 'moment_secondary_body', moment_secondary);
if options.geometry == "planar"
    differential = rotation_secondary*force_secondary - ...
        force_primary_inertial*(params.mu_Dimorphos/params.mu_Didymos);
    radiation.planar_force = [dot(differential, radial); dot(differential, tangential)];
end
radiation.diagnostics = struct( ...
    'maximum_shadowed_faces', nnz(shadow_primary)+nnz(shadow_secondary)+ ...
        double(eclipse < 1)*numel(secondary.area), ...
    'maximum_visible_faces', numel(visible_primary)+numel(visible_secondary), ...
    'maximum_mutual_flux_W_m2', max([mutual_primary;mutual_secondary]), ...
    'maximum_recoil_N', norm(force_secondary), ...
    'maximum_rough_recoil_difference_N', rough_difference);
end

function mu = local_incidence(body, direction, external_shadow)
mu = max(0, body.normal*direction);
mu(external_shadow) = 0;
faces = find(body.shadow.occ_rev_map > 0 & mu > 0);
bin = cobird.geometry.compute_sun_direction_bin(direction, body.shadow);
blocked = body.shadow.shadow_table(double(body.shadow.occ_rev_map(faces)), bin);
mu(faces(blocked)) = 0;
end

function source = local_sources(body, temperature, mu, direction, params, distance_au, model)
rough = model.options.surface_mode == "rough";
source = struct('mu', mu, 'emission', params.epsilon*params.sigma*temperature.smooth(:, 1).^4, ...
    'emission_rough', [], 'direct_rough', [], 'reflected_rough', []);
source.reflected = cobird.thermal.compute_reflected_light_flux(params.A, params.F_s/distance_au^2*mu);
if rough
    source.emission_rough = cobird.thermal.rough_surface_flux( ...
        temperature.rough, model.crater, numel(body.area), params);
    source.direct_rough = cobird.thermal.compute_rough_direct_solar_flux( ...
        body.normal, body.x_local, body.y_local, direction, params, distance_au, model.crater, mu);
    source.reflected_rough = cobird.thermal.compute_reflected_light_flux(params.A, source.direct_rough);
end
[source.extra, source.extra_rough] = cobird.thermal.compute_surface_exchange( ...
    body.rough_model, body.exchange, source.emission, source.emission_rough, ...
    source.reflected, source.reflected_rough, params, rough, model.options.terrain_heating);
end

function [force, moment, difference] = local_recoil(body, temperature, source, params, rough)
[smooth_force, smooth_moment] = cobird.thermal.compute_radiative_recoil( ...
    body.normal, body.centroid, body.area, temperature.smooth(:, 1), params, ...
    source.reflected, body.occultation);
force = smooth_force;
moment = smooth_moment;
difference = 0;
if rough
    [rough_force, rough_moment] = cobird.thermal.compute_directional_rough_recoil( ...
        body.rough_model.recoil, source.emission_rough+source.reflected_rough);
    force = (1-params.roughness_f_R)*smooth_force+rough_force;
    moment = (1-params.roughness_f_R)*smooth_moment+rough_moment;
    difference = norm(force-smooth_force);
end
end
