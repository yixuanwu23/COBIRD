function operator = finalize_rough_recoil_operator( ...
        operator, centroid, area, fraction, light_speed)
%FINALIZE_ROUGH_RECOIL_OPERATOR Contract fixed geometry into a force map.

assert(strcmp(operator.schema, 'COBIRD_ROUGH_RECOIL_INTEGRAL'), ...
    'cobird:geometry:finalize_rough_recoil_operator:schema', ...
    'The angular recoil integral is incompatible.');
parent_count = numel(area);
crater_count = size(operator.local_coefficients, 1);
if isscalar(fraction)
    fraction = repmat(fraction, parent_count, 1);
end
fraction = fraction(:);
assert(isequal(size(centroid), [parent_count, 3]) && ...
    numel(fraction) == parent_count && ...
    all(fraction >= 0 & fraction <= 1) && ...
    isfinite(light_speed) && light_speed > 0, ...
    'cobird:geometry:finalize_rough_recoil_operator:input', ...
    'Centroids, areas, rough fractions, or light speed are invalid.');

coefficients = zeros(crater_count, parent_count, 3);
frame = cat(3, operator.x_axis, operator.y_axis, operator.normal);
for component = 1:3
    for axis = 1:3
        coefficients(:, :, component) = coefficients(:, :, component) + ...
            operator.local_coefficients(:, :, axis) .* frame(:, component, axis)';
    end
    coefficients(:, :, component) = coefficients(:, :, component) .* ...
        (-area(:) .* fraction / light_speed)';
end

force_map = reshape(coefficients, crater_count * parent_count, 3)';
lever = repelem(centroid, crater_count, 1);
torque_map = cross(lever, force_map', 2)';
operator.net_map = [force_map; torque_map];
operator.force_coefficients = coefficients;
operator.centroid = centroid;
operator.parent_count = parent_count;
operator.crater_count = crater_count;
operator.fraction = fraction;
operator.schema = 'COBIRD_ROUGH_RECOIL_OPERATOR';
operator = rmfield(operator, ...
    {'local_coefficients', 'normal', 'x_axis', 'y_axis'});
end
