function [indices, coefficients, active] = crater_direction_stencil(table, direction)
%CRATER_DIRECTION_STENCIL Return bilinear table indices and coefficients.

assert(size(direction, 2) == 3 && all(isfinite(direction), 'all') && ...
    all(abs(vecnorm(direction, 2, 2) - 1) < 1e-8), ...
    'cobird:geometry:crater_direction_stencil:direction', ...
    'Directions must be finite unit vectors.');
assert(strcmp(table.schema, 'COBIRD_CRATER_DIRECTION_TABLE'), ...
    'cobird:geometry:crater_direction_stencil:schema', ...
    'A COBIRD equal-angle crater-direction table is required.');

direction_count = size(direction, 1);
indices = ones(4, direction_count);
coefficients = zeros(4, direction_count);
active = direction(:, 3) > 0;
direction = direction(active, :);

azimuth = mod(atan2(direction(:, 2), direction(:, 1)), 2 * pi) * ...
    table.azimuth_count / (2 * pi);
% A rounded full turn belongs to the first azimuth bin.
azimuth = mod(azimuth, table.azimuth_count);
azimuth_0 = floor(azimuth);
azimuth_fraction = azimuth - azimuth_0;
azimuth_1 = mod(azimuth_0 + 1, table.azimuth_count);
elevation = asin(max(0, min(1, direction(:, 3)))) / table.elevation_step;
elevation = min(table.elevation_count - 1, elevation);
elevation_0 = min(table.elevation_count - 2, floor(elevation));
elevation_fraction = elevation - elevation_0;

column_00 = 1 + azimuth_0 + elevation_0 * table.azimuth_count;
column_10 = 1 + azimuth_1 + elevation_0 * table.azimuth_count;
indices(:, active) = [column_00, column_10, ...
    column_00 + table.azimuth_count, ...
    column_10 + table.azimuth_count]';
coefficients(:, active) = [ ...
    (1 - azimuth_fraction) .* (1 - elevation_fraction), ...
    azimuth_fraction .* (1 - elevation_fraction), ...
    (1 - azimuth_fraction) .* elevation_fraction, ...
    azimuth_fraction .* elevation_fraction]';
end
