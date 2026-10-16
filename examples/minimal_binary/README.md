# Minimal directional-roughness example

This data-free example uses two generated polyhedra with roughness on both
bodies. The parent meshes are small; crater shadowing and recoil use the core's
1-degree angular resolution.

With a supported C++ compiler configured, build the two platform-specific
radiation kernels, then run from the repository root:

```matlab
startup_cobird
cobird.thermal.build_radiation_mex
run(fullfile("examples", "minimal_binary", "run_example.m"))
```

No external shape model or SPICE kernel is required. Temporary crater data
are removed when the function exits; no result file is retained. The default
run advances 48 steps of 300 seconds. It prints checks and leaves the returned
`result` in the workspace, with diagnostics in `result.checks` and sampled
outputs in `result.solution`.

The example function accepts surface and orbit options:

```matlab
addpath(fullfile("examples", "minimal_binary"))
result = run_minimal_binary(surface_mode="rough", coupled=true);
```

Set `surface_mode="smooth"` to use smooth facets; this mode needs no MEX kernels.
With `coupled=false`, a prescribed circular orbit supplies the geometry.
With `coupled=true`, the orbit and libration evolve under the computed recoil.
The primary attitude is prescribed in both modes.
The convex parent meshes need no terrain self-shadowing or self-heating.
Every depth node starts at 220 K. The example uses explicit thermal steps;
the implicit solver is available through `cobird.prepare_binary`.
