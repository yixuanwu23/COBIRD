function report = check_dependencies(options)
%CHECK_DEPENDENCIES Report COBIRD, MICE, compiler, and MEX status.

arguments
    options.require_mice (1,1) logical = false
    options.require_thermal_mex (1,1) logical = false
    options.require_radiation_mex (1,1) logical = false
end

setup_root = fileparts(mfilename('fullpath'));
project_root = fileparts(setup_root);
addpath(fullfile(project_root, 'src'));

report = struct();
report.matlab_version = version;
report.matlab_release = version('-release');
report.mice_root = cobird.runtime.setup_project_paths( ...
    project_root, require_mice=options.require_mice);
report.mice_version = '';
if ~isempty(report.mice_root)
    report.mice_version = cspice_tkvrsn('toolkit');
end

try
    compiler = mex.getCompilerConfigurations('C++', 'Selected');
    if isempty(compiler)
        report.cpp_compiler = '';
    else
        report.cpp_compiler = compiler.Name;
    end
catch
    report.cpp_compiler = '';
end

required_mex = { ...
    'cobird.thermal.calculate_temperature_tiled_threads_mex', ...
    'cobird.thermal.calculate_temperature_fixed20_threads_mex', ...
    'cobird.thermal.calculate_temperature_nonuniform_threads_mex'};
report.thermal_mex_available = all(cellfun(@(name) ~isempty(which(name)), required_mex));
if options.require_thermal_mex && ~report.thermal_mex_available
    error('cobird:setup:missingThermalMex', ...
        'Thermal MEX solvers are missing. Run cobird.thermal.build_mex_solvers.');
end

radiation_mex = { ...
    'cobird.geometry.interpolate_crater_weights_mex', ...
    'cobird.thermal.contract_mutual_direction_mex'};
report.radiation_mex_available = ...
    all(cellfun(@(name) ~isempty(which(name)), radiation_mex));
if options.require_radiation_mex && ~report.radiation_mex_available
    error('cobird:setup:missingRadiationMex', ...
        ['Directional-radiation MEX kernels are missing. ' ...
         'Run cobird.thermal.build_radiation_mex.']);
end

fprintf('MATLAB:      %s (%s)\n', report.matlab_version, report.matlab_release);
if isempty(report.mice_root)
    fprintf('MICE:        not configured\n');
    fprintf('Download:    https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html\n');
else
    fprintf('MICE:        %s\n', report.mice_version);
end
if isempty(report.cpp_compiler)
    fprintf('C++ compiler: not configured\n');
else
    fprintf('C++ compiler: %s\n', report.cpp_compiler);
end
fprintf('Thermal MEX: %s\n', local_status(report.thermal_mex_available));
fprintf('Radiation MEX: %s\n', local_status(report.radiation_mex_available));
end

function status = local_status(value)
if value
    status = 'available';
else
    status = 'not built';
end
end
