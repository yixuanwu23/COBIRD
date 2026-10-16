function weights = evaluate_crater_direction_table(table, direction)
%EVALUATE_CRATER_DIRECTION_TABLE Interpolate nonnegative crater weights.
% Directions on or below the aperture horizon transmit no power.

[indices, coefficients, active] = ...
    cobird.geometry.crater_direction_stencil(table, direction);
if ~any(active)
    weights = zeros(size(table.weights, 1), size(direction, 1));
    return
end
weights = cobird.geometry.interpolate_crater_weights_mex( ...
    table.weights, indices, coefficients);
end
