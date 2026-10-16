function flux = rough_surface_flux(temperature, crater, parent_count, params)
%ROUGH_SURFACE_FLUX Return the exitance of every crater facet.

crater_count = size(crater.normal, 1);
assert(~isempty(temperature) && ...
    size(temperature, 1) == crater_count * parent_count, ...
    'cobird:thermal:rough_surface_flux:state', ...
    'The crater temperature state is missing or incompatible.');
surface_temperature = reshape(temperature(:, 1), crater_count, parent_count);
flux = params.epsilon * params.sigma * surface_temperature.^4;
end
