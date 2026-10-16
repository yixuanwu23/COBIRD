function [config, provenance] = resolve_run_configuration(defaults, prefix)
%RESOLVE_RUN_CONFIGURATION Apply named environment overrides to defaults.
% Environment names follow field paths. For prefix='RUN_':
%   forces.tidal   -> RUN_FORCES_TIDAL
%   duration.days  -> RUN_DURATION_DAYS
% JSON overrides are intentionally unsupported.
arguments
    defaults (1, 1) struct
    prefix (1, :) char = 'RUN_'
end

config = defaults;
overrides = struct();
[config, overrides] = local_apply(config, defaults, prefix, '', overrides);

provenance = struct();
provenance.precedence = 'code_defaults < named_environment';
provenance.code_defaults = defaults;
provenance.named_environment_overrides = overrides;
provenance.effective_config = config;
end

function [out, overrides] = local_apply(out, defaults, prefix, path, overrides)
names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    value = defaults.(name);
    if isempty(path)
        field_path = name;
    else
        field_path = [path, '.', name];
    end
    if isstruct(value) && isscalar(value)
        [out.(name), overrides] = local_apply(out.(name), value, prefix, field_path, overrides);
        continue
    end
    env_name = upper([prefix, strrep(field_path, '.', '_')]);
    raw = strtrim(getenv(env_name));
    if isempty(raw)
        continue;
    end
    parsed = local_parse(raw, value, env_name);
    out.(name) = parsed;
    overrides.(matlab.lang.makeValidName(field_path)) = ...
        struct('environment', env_name, 'value', parsed);
end
end

function value = local_parse(raw, prototype, env_name)
if islogical(prototype) && isscalar(prototype)
    if any(strcmpi(raw, {'1', 'true', 'yes', 'on'}))
        value = true;
    elseif any(strcmpi(raw, {'0', 'false', 'no', 'off'}))
        value = false;
    else
        error('cobird:runtime:resolve_run_configuration:badLogical', ...
            '%s must be true/false, yes/no, on/off, or 1/0.', env_name);
    end
elseif isnumeric(prototype) && (isscalar(prototype) || isempty(prototype))
    value = str2double(raw);
    assert(isfinite(value), 'cobird:runtime:resolve_run_configuration:badNumber', ...
        '%s must be a finite number.', env_name);
elseif ischar(prototype) || (isstring(prototype) && isscalar(prototype))
    value = raw;
else
    error('cobird:runtime:resolve_run_configuration:unsupportedType', ...
        '%s cannot override configuration values of type %s.', ...
        env_name, class(prototype));
end
end
