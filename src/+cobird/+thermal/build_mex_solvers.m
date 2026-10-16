function build_mex_solvers()
%BUILD_MEX_SOLVERS Compile thermal-column and directional-radiation MEX solvers.

package_dir = fileparts(mfilename('fullpath'));
source_files = { ...
    'calculate_temperature_tiled_threads_mex.cpp', ...
    'calculate_temperature_fixed20_threads_mex.cpp', ...
    'calculate_temperature_nonuniform_threads_mex.cpp'};

if ispc
    compiler_flags = 'COMPFLAGS="$COMPFLAGS /utf-8 /O2"';
    linker_flags = {};
elseif isunix
    compiler_flags = 'CXXFLAGS="$CXXFLAGS -O3 -std=c++11 -pthread"';
    linker_flags = {'LDFLAGS="$LDFLAGS -pthread"'};
else
    error('cobird:thermal:build_mex_solvers:unsupportedPlatform', ...
        'The thermal-column MEX build script does not support this platform.');
end

for source_index = 1:numel(source_files)
    source_file = fullfile(package_dir, source_files{source_index});
    assert(isfile(source_file), 'cobird:thermal:build_mex_solvers:missingSource', ...
        'MEX source file not found: %s', source_file);
    mex('-R2018a', compiler_flags, linker_flags{:}, ...
        '-outdir', package_dir, source_file);
end
cobird.thermal.build_radiation_mex();
fprintf('Thermal and radiation MEX solvers built successfully for .%s.\n', mexext);
end
