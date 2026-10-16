function iv = integrate_terrain_horizon(v1, v2, v3)
%INTEGRATE_TERRAIN_HORIZON Integrate recoil above the local terrain horizon.
% Each input is an M-by-3 array of triangle vertices in the source facet's
% frame: origin at its centroid, z along its normal. The highest triangle
% intersection defines the horizon; directions below it are treated as blocked.
if isempty(v1)
    iv = [0, 0, 2*pi];
    return
end

x = [v1(:, 1), v2(:, 1), v3(:, 1)];
y = [v1(:, 2), v2(:, 2), v3(:, 2)];
z = [v1(:, 3), v2(:, 3), v3(:, 3)];
tolerance = 32 * eps(max(abs([x, y, z]), [], 2));

% A triangle strictly above the origin closes the filled horizon at all azimuths.
scale = max(abs([x, y]), [], 2);
rows = find(scale > 0);
px = x(rows, :) ./ scale(rows);
py = y(rows, :) ./ scale(rows);
next = [2, 3, 1];
signed_area = px .* py(:, next) - py .* px(:, next);
inside = all(signed_area > 64*eps, 2) | all(signed_area < -64*eps, 2);
rows = rows(inside);
signed_area = signed_area(inside, :);
weights = signed_area(:, [2, 3, 1]) ./ sum(signed_area, 2);
axis_z = sum(weights .* z(rows, :), 2);
if any(axis_z > tolerance(rows))
    iv = [0, 0, 0];
    return
end

% Split quadrature at angular support boundaries, including z=0 crossings.
angles = atan2(y(:), x(:));
for edge = 1:3
    next = mod(edge, 3) + 1;
    crosses = (z(:, edge) < 0 & z(:, next) > 0) | ...
        (z(:, edge) > 0 & z(:, next) < 0);
    t = z(crosses, edge) ./ (z(crosses, edge) - z(crosses, next));
    px = x(crosses, edge) + t .* (x(crosses, next) - x(crosses, edge));
    py = y(crosses, edge) + t .* (y(crosses, next) - y(crosses, edge));
    angles = [angles; atan2(py, px)]; %#ok<AGROW>
end
waypoints = unique(mod(angles, 2*pi));
waypoints = waypoints(waypoints > 0 & waypoints < 2*pi);
iv = [0, 0, 2*pi];
for component = 1:3
    iv(component) = iv(component) + integral(@(phi) local_integrand( ...
        phi, component, x, y, z, tolerance), 0, 2*pi, ...
        'Waypoints', waypoints, 'AbsTol', 1e-10, 'RelTol', 1e-8);
end
end

function value = local_integrand(phi, component, x, y, z, tolerance)
value = zeros(size(phi));
% Bound temporary triangle-by-angle arrays for large candidate sets.
batch_size = max(1, floor(32768 / size(x, 1)));
for first = 1:batch_size:numel(phi)
    indices = first:min(first + batch_size - 1, numel(phi));
    angles = reshape(phi(indices), 1, []);
    height = local_horizon(angles, x, y, z, tolerance);
    if component == 1
        value(indices) = (cos(height).^3 - 1) .* cos(angles);
    elseif component == 2
        value(indices) = (cos(height).^3 - 1) .* sin(angles);
    else
        value(indices) = -sin(height).^3;
    end
end
end

function height = local_horizon(phi, x, y, z, tolerance)
count = size(x, 1);
rho = cell(1, 3);
side = cell(1, 3);
points_rho = nan(count, numel(phi), 6);
points_z = nan(size(points_rho));
vertical = false(count, numel(phi));
for vertex = 1:3
    rho{vertex} = x(:, vertex) * cos(phi) + y(:, vertex) * sin(phi);
    side{vertex} = y(:, vertex) * cos(phi) - x(:, vertex) * sin(phi);
    on_plane = abs(side{vertex}) <= tolerance;
    r = rho{vertex};
    r(~on_plane) = NaN;
    points_rho(:, :, vertex) = r;
    points_z(:, :, vertex) = repmat(z(:, vertex), 1, numel(phi));
end

for edge = 1:3
    next = mod(edge, 3) + 1;
    a = side{edge};
    b = side{next};
    crosses = (a < 0 & b > 0) | (a > 0 & b < 0);
    t = nan(size(a));
    t(crosses) = a(crosses) ./ (a(crosses) - b(crosses));
    points_rho(:, :, edge + 3) = rho{edge} + t .* (rho{next} - rho{edge});
    points_z(:, :, edge + 3) = z(:, edge) + t .* (z(:, next) - z(:, edge));

    % A coplanar edge can cross the upward axis between its endpoints.
    on_plane = abs(a) <= tolerance & abs(b) <= tolerance;
    crosses_axis = on_plane & ((rho{edge} <= tolerance & rho{next} > tolerance) | ...
        (rho{next} <= tolerance & rho{edge} > tolerance));
    t(:) = NaN;
    t(crosses_axis) = rho{edge}(crosses_axis) ./ ...
        (rho{edge}(crosses_axis) - rho{next}(crosses_axis));
    axis_z = z(:, edge) + t .* (z(:, next) - z(:, edge));
    vertical = vertical | (crosses_axis & axis_z > tolerance);
end

elevation = atan2(points_z, points_rho);
elevation(~(points_rho > tolerance & points_z > tolerance)) = 0;
height = max(elevation, [], 3);
[r_min, i_min] = min(points_rho, [], 3, 'omitnan');
[r_max, i_max] = max(points_rho, [], 3, 'omitnan');
offset = reshape(1:numel(r_min), size(r_min));
z_min = points_z(offset + (i_min - 1) * numel(r_min));
z_max = points_z(offset + (i_max - 1) * numel(r_min));
crosses_axis = r_min <= tolerance & r_max > tolerance;
axis_z = zeros(size(r_min));
axis_z(crosses_axis) = z_min(crosses_axis) - r_min(crosses_axis) .* ...
    (z_max(crosses_axis) - z_min(crosses_axis)) ./ ...
    (r_max(crosses_axis) - r_min(crosses_axis));
height(vertical | (crosses_axis & axis_z > tolerance)) = pi/2;
height = max(height, [], 1);
end
