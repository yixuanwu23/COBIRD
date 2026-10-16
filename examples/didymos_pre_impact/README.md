# Didymos pre-impact example

This is the primary COBIRD example. It loads the official NASA PDS Didymos
v003 and final pre-impact Dimorphos v004 meshes without modifying either
source file. The example exercises eclipses, mutual infrared heating,
smooth-surface thermal evolution, radiative recoil, and planar
orbit-libration coupling.

The model is assembled with `cobird.prepare_binary` and evolved with
`cobird.run_coupled`.

Run from the repository root:

```matlab
startup_cobird
run(fullfile("examples", "didymos_pre_impact", "run_example.m"))
```

The default run advances 24 steps of 300 seconds. Its checks cover source
integrity, finite outputs, and activation of the main calculations.
`result.checks` contains these checks, and `result.solution`
contains sampled states, temperatures, forces, and torques. For example,
`result.solution.final_temperature.secondary.smooth(:,1)` contains the final
secondary surface temperatures in kelvin. No result file is written.
MICE, SPICE kernels, and compiled MEX solvers are not required.

To set the step count or time step, call the example function:

```matlab
addpath(fullfile("examples", "didymos_pre_impact"))
result = run_didymos_pre_impact(step_count=24, time_step_s=300);
```

## Model assumptions

- The OBJ and XML files are the unmodified PDS products under
  `data/shape_models/dart_pre_impact/`.
- Vertex coordinates are converted from kilometres to metres in memory.
  Both meshes remain fixed throughout the run.
- The pre-impact mutual period and semimajor axis set the total mass. The
  bodies have equal bulk density, so their masses are proportional to mesh volume.
- Gravity coefficients and secondary inertia use ellipsoid approximations
  whose semiaxes are half the coordinate extents of the meshes.
- The Sun is fixed at 1 au; the primary rotates at a prescribed 2.2600-hour period.
- Every thermal depth node starts at 200 K. The two-hour run shows the initial
  thermal response to illumination.
- Tides and the empirical DART drag term are disabled.

Macroscopic terrain self-shadowing and terrain radiative exchange are disabled.
The calculation uses the full PDS meshes with smooth surfaces and explicit
thermal steps. Directional roughness is demonstrated in
[minimal_binary](../minimal_binary/README.md).
