function [config, provenance, primary, secondary] = resolve_separate_tidal_parameters(config, provenance)
%RESOLVE_SEPARATE_TIDAL_PARAMETERS Resolve shared and body-specific k2/Q.
% Explicit body values override the shared value; named environment values
% override both. No environment variables are changed by this function.

primary = config.tidal.k2_over_Q;
secondary = primary;
names = {'k2_over_Q_primary', 'k2_over_Q_secondary'};
for index = 1:numel(names)
    field = names{index};
    value = config.tidal.k2_over_Q;
    if isfield(config.tidal, field)
        value = config.tidal.(field);
    end
    environment = upper(['RUN_TIDAL_', field]);
    raw = strtrim(getenv(environment));
    if ~isempty(raw)
        value = str2double(raw);
        provenance.named_environment_overrides.(['tidal_', field]) = ...
            struct('environment', environment, 'value', value);
    end
    assert(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0, ...
    'cobird:runtime:resolve_separate_tidal_parameters:invalidValue', ...
        '%s must be positive and finite.', field);
    config.tidal.(field) = value;
    if index == 1
        primary = value;
    else
        secondary = value;
    end
end
provenance.effective_config = config;
end
