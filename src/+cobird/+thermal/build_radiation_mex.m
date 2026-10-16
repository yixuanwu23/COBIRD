function build_radiation_mex()
%BUILD_RADIATION_MEX Compile directional interpolation and exchange kernels.

thermal_dir = fileparts(mfilename('fullpath'));
geometry_dir = fullfile(fileparts(thermal_dir), '+geometry');
if ispc
    compiler_flags = 'COMPFLAGS="$COMPFLAGS /utf-8 /O2 /fp:strict"';
elseif isunix
    compiler_flags = ...
        'CXXFLAGS="$CXXFLAGS -O3 -std=c++11 -ffp-contract=off"';
else
    error('cobird:thermal:build_radiation_mex:unsupportedPlatform', ...
        'The directional-radiation MEX build is not supported on this platform.');
end
mex('-R2018a', compiler_flags, '-outdir', geometry_dir, ...
    fullfile(geometry_dir, 'interpolate_crater_weights_mex.cpp'));
mex('-R2018a', compiler_flags, '-outdir', thermal_dir, ...
    fullfile(thermal_dir, 'contract_mutual_direction_mex.cpp'));
fprintf('Directional-radiation MEX kernels built successfully for .%s.\n', mexext);
end
