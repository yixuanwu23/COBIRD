function Q_direct = compute_rough_direct_solar_flux( ...
    normal_parent, x_parent, y_parent, dir_sun_body, params, ...
    d_sun_AU, crater, mu_parent)
%COMPUTE_ROUGH_DIRECT_SOLAR_FLUX Direct solar flux on explicit crater facets.

parent_count = size(normal_parent, 1);
crater_count = size(crater.normal, 1);
assert(numel(mu_parent) == parent_count, ...
    'cobird:thermal:compute_rough_direct_solar_flux:parentSize', ...
    'mu_parent size does not match the parent mesh.');

sun_parent = [x_parent * dir_sun_body, ...
              y_parent * dir_sun_body, ...
              normal_parent * dir_sun_body];
mu_crater = max(0, crater.normal * sun_parent');

parent_cosine = max(0, sun_parent(:, 3));
illumination_scale = zeros(parent_count, 1);
active_parent = parent_cosine > 0 & mu_parent(:) > 0;
illumination_scale(active_parent) = ...
    mu_parent(active_parent) ./ parent_cosine(active_parent);
mu_crater = mu_crater .* illumination_scale';

if any(active_parent)
    sx = sun_parent(active_parent, 1);
    sy = sun_parent(active_parent, 2);
    sz = sun_parent(active_parent, 3);
    az = atan2(sy, sx);
    el = asin(max(-1, min(1, sz)));
    az_bin = max(1, min(crater.N_az_shadow, ...
        floor((az + pi) / crater.daz) + 1));
    el_bin = max(1, min(crater.N_el_shadow, ...
        floor((el + pi/2) / crater.del) + 1));
    bin = (el_bin - 1) * crater.N_az_shadow + az_bin;

    mu_active = mu_crater(:, active_parent);
    mu_active(crater.shadow_table(:, bin)) = 0;
    mu_crater(:, active_parent) = mu_active;
end

assert(isequal(size(mu_crater), [crater_count, parent_count]));
Q_direct = params.F_s / d_sun_AU^2 .* mu_crater;
end
