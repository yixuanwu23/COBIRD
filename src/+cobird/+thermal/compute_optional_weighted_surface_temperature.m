function mean_temperature = compute_optional_weighted_surface_temperature( ...
        macro_temperature, rough_temperature, crater, params, use_explicit_roughness)
%COMPUTE_OPTIONAL_WEIGHTED_SURFACE_TEMPERATURE Blend available rough temperatures.

macro_surface = macro_temperature(:, 1);
if ~use_explicit_roughness || isempty(rough_temperature)
    mean_temperature = macro_surface;
    return;
end

crater_count = size(crater.normal, 1);
parent_count = numel(macro_surface);
rough_surface = reshape(rough_temperature(:, 1), crater_count, parent_count);
rough_mean = (crater.area(:)' * rough_surface)' / sum(crater.area);
rough_fraction = local_expand_fraction(params.roughness_f_R, parent_count);
mean_temperature = (1 - rough_fraction) .* macro_surface + ...
    rough_fraction .* rough_mean;
end

function fraction = local_expand_fraction(fraction, count)
if isscalar(fraction)
    fraction = repmat(fraction, count, 1);
else
    fraction = fraction(:);
end
end
