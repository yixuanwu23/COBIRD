function weights = crater_directional_weights(crater, direction)
%CRATER_DIRECTIONAL_WEIGHTS Project visible crater walls onto the aperture.
% The result is normalized by the actual mesh aperture and is used for
% both directional emission and reciprocal reception.

validateattributes(direction, {'double'}, {'2d', 'ncols', 3, 'finite'});
assert(all(abs(vecnorm(direction, 2, 2) - 1) < 1e-10), ...
    'cobird:geometry:crater_directional_weights:notUnit', ...
    'Directions must be unit vectors.');

vertices = crater.vertices;
rim = vertices(abs(vertices(:, 3)) < 1e-10, 1:2);
assert(size(rim, 1) >= 3, ...
    'cobird:geometry:crater_directional_weights:noRim', ...
    'The crater mesh has no aperture rim at z = 0.');
[~, order] = sort(atan2(rim(:, 2), rim(:, 1)));
rim = rim(order, :);
edge = rim([2:end, 1], :) - rim;
offset = edge(:, 2) .* rim(:, 1) - edge(:, 1) .* rim(:, 2);

crater_count = size(crater.faces, 1);
weights = zeros(crater_count, size(direction, 1));
for direction_index = 1:size(direction, 1)
    ray = direction(direction_index, :);
    if ray(3) <= 0
        continue
    end
    planes = [-edge(:, 2) * ray(3), edge(:, 1) * ray(3), ...
        edge(:, 2) * ray(1) - edge(:, 1) * ray(2)];
    constants = offset * ray(3);
    cosine = crater.normal * ray';
    for crater_index = find(cosine > 0)'
        polygon = vertices(crater.faces(crater_index, :), :);
        for rim_index = 1:size(rim, 1)
            polygon = local_clip(polygon, planes(rim_index, :), constants(rim_index));
            if size(polygon, 1) < 3
                break
            end
        end
        if size(polygon, 1) < 3
            continue
        end
        area_vector = sum(cross(polygon, polygon([2:end, 1], :), 2), 1) / 2;
        visible_area = abs(dot(area_vector, crater.normal(crater_index, :)));
        weights(crater_index, direction_index) = visible_area * ...
            cosine(crater_index) / (crater.opening_area * ray(3));
    end
end
end

function clipped = local_clip(polygon, normal, offset)
value = polygon * normal' + offset;
inside = value >= 0;
if all(inside)
    clipped = polygon;
    return
elseif ~any(inside)
    clipped = zeros(0, 3);
    return
end

clipped = zeros(size(polygon, 1) + 1, 3);
count = 0;
previous = size(polygon, 1);
for current = 1:size(polygon, 1)
    if inside(previous) ~= inside(current)
        fraction = value(previous) / (value(previous) - value(current));
        count = count + 1;
        clipped(count, :) = polygon(previous, :) + fraction * ...
            (polygon(current, :) - polygon(previous, :));
    end
    if inside(current)
        count = count + 1;
        clipped(count, :) = polygon(current, :);
    end
    previous = current;
end
clipped = clipped(1:count, :);
end
