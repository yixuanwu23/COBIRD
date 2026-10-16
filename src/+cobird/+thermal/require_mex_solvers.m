function require_mex_solvers(options)
%REQUIRE_MEX_SOLVERS Check the compiled solvers needed by a workflow.

arguments
    options.thermal (1, 1) logical = true
    options.radiation (1, 1) logical = true
end

required = { ...
    'cobird.thermal.calculate_temperature_tiled_threads_mex', ...
    'cobird.thermal.calculate_temperature_fixed20_threads_mex', ...
    'cobird.thermal.calculate_temperature_nonuniform_threads_mex'};
if ~options.thermal
    required = {};
end
if options.radiation
    required = [required, { ...
        'cobird.geometry.interpolate_crater_weights_mex', ...
        'cobird.thermal.contract_mutual_direction_mex'}];
end
missing = required(cellfun(@(name) isempty(which(name)), required));
if ~isempty(missing)
    error('cobird:thermal:require_mex_solvers:missingSolver', ...
        ['Required MEX solvers are missing for this platform (%s). ' ...
         'Run cobird.thermal.build_mex_solvers first.'], strjoin(missing, ', '));
end
end
