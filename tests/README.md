# Tests

The full suite requires MATLAB R2024a or later and a supported C++ compiler.
Run from the repository root; the build command compiles both thermal and
directional-radiation MEX sources for the current platform:

```matlab
startup_cobird
cobird.thermal.build_mex_solvers
results = runtests("tests");
assertSuccess(results)
```

The core tests generate geometry in memory or in a temporary directory. The
mutual-visibility tests check positive-cosine retention, forward and reverse
exchange, body-frame shadow queries, inverse-square scaling, and invalid
geometry inputs. A separate data-integrity test verifies the bundled PDS
OBJ/XML hashes. The suite does not require MICE or SPICE kernels. The synthetic
directional-roughness test uses MEX kernels rebuilt locally from the included
C++ sources; compiled binaries are not distributed.

Short physical-driver checks cover smooth/rough thermal and coupled runs,
3D state propagation, output cadence, independent tidal parameters, and
checkpoint continuation with the actual crater-file fingerprint.
