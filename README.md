# COBIRD

**COBIRD** (**CO**upled **BI**nary **R**adiative **D**ynamics) is MATLAB research
software that evolves temperature profiles on both bodies of a binary asteroid,
their relative orbit, and the secondary's attitude. It includes eclipses, mutual
infrared heating, hemispherical-crater roughness, thermal and reflected-light
recoil, and planar or three-dimensional dynamics. The primary attitude is
prescribed by an environment callback.

This repository accompanies the IAC study *Post-Impact Thermophysical-Dynamical
Evolution of the Didymos-Dimorphos System under Coupled Yarkovsky-YORP Effects*.

## Scope

The repository contains the reusable numerical core, the official NASA PDS
Didymos v003 and final pre-impact Dimorphos v004 shape products, and two runnable
examples.

## Requirements

MATLAB R2024a or later is required. The PDS quick example needs no compiler,
MICE, or SPICE kernels. Directional roughness, implicit thermal solvers, and
the full test suite require a supported C++ compiler and locally built MEX
files. Build both thermal and radiation solvers with
`cobird.thermal.build_mex_solvers` after `startup_cobird`.

[NAIF MICE](external/mice/README.md) and user-supplied
[SPICE kernels](data/kernels/README.md) are needed only for SPICE callbacks.
After startup, `check_dependencies` reports the configured dependencies.

## Quick start

Start MATLAB in the repository root and run:

```matlab
startup_cobird
run(fullfile("examples", "didymos_pre_impact", "run_example.m"))
```

The primary example reads the original PDS meshes, verifies their SHA-256
values, and advances 24 coupled steps of 300 seconds with smooth surfaces.
It prints checks for source integrity, finite outputs, and activation of the
example's main calculations, followed by a temperature range. It leaves `result`
in the workspace. Inspect `result.checks`, `result.solution.time_s`, and
`result.solution.final_temperature`; no result file is written. Kilometre
coordinates are converted to metres only in memory. See the
[example README](examples/didymos_pre_impact/README.md) for its assumptions.

The data-free `examples/minimal_binary/` case exercises directional roughness
on both bodies of a small synthetic binary. Build the directional
radiation kernels before running it:

```matlab
cobird.thermal.build_radiation_mex
run(fullfile("examples", "minimal_binary", "run_example.m"))
```

Both examples use `cobird.prepare_binary` and the time-stepping drivers.
See [Tests](tests/README.md) for the full suite and build command.

## Layout

```text
COBIRD/
|-- src/+cobird/              Namespaced MATLAB and C++ source
|-- examples/                PDS pre-impact and synthetic examples
|-- tests/                   Unit, integration, and data-integrity tests
|-- setup/                   Dependency checks
|-- data/kernels/             Kernel manifest and meta-kernel templates
|-- data/shape_models/        Original PDS pre-impact OBJ/XML products
|-- external/mice/            Installation instructions only
|-- LICENSE                  MIT license for COBIRD
|-- LICENSES/                Third-party license texts
|-- README.md
|-- CITATION.cff
`-- THIRD_PARTY_NOTICES.md
```

## Calculation interfaces

`cobird.prepare_binary(primary, secondary, params, ...)` assembles heating,
recoil, and dynamics for two fixed meshes. The example functions define complete
inputs; `help cobird.prepare_binary` describes the input fields and conventions.
The main choices are:

- `surface_mode="smooth"` or `"rough"`; rough mode requires `crater_file`.
- `geometry="planar"` (6 states): `[r; rdot; theta; thetadot; phi; phidot]`, with
  secondary libration angle `phi` relative to the line of centers.
- `geometry="3d"` (13 states): relative position and velocity, scalar-first
  secondary-body-to-inertial quaternion, and secondary body-frame angular velocity.
- `thermal_solver="implicit"` (default, compiled solvers) or `"explicit"`
  (uniform depth grids and a suitably small time step).
- `terrain_heating=true` (default) enables radiative exchange within each body
  and requires its view-factor cache from `build_view_factor_cache`.
  Geometry and caches must use metres; see the
  [PDS data README](data/shape_models/dart_pre_impact/README.md).
- `mutual_heating=true` (default) enables infrared heating of the secondary by
  the primary. `reverse_heating=true` adds heating of the primary by the secondary;
  its default is `false`.
- Each body can supply its own `dz` and `Nz`; otherwise `params.dz` and
  `params.Nz` apply to both. `dz` may be a spacing or a depth-node vector.

The relative position points from primary to secondary. In planar mode,
`phi=0` points the secondary's +x axis toward the primary and its +z axis along
the orbit normal. The columns of `orbit_basis` are the radial direction at
`theta=0`, the tangential direction, and the orbit normal. In 3D mode, `params.I`
is the secondary body-frame inertia tensor divided by its mass, in m².
`forces.secondary_gravity` controls its quadrupole orbital acceleration and
gravity-gradient torque together and defaults to `true`.

`environment(t)` supplies primary-relative `sun_position_m` and
`primary_to_inertial` in the same inertial frame chosen for the state and
`orbit_basis`. Enabled 3D
perturbers also need primary-relative `EMB_position_m`, `Venus_position_m`,
or `Jupiter_position_m`. Time is in seconds, positions in metres, angles in
radians, and rotations map body coordinates to the inertial frame.
Initial and prescribed quaternions are normalized. Planar radius must be
positive, and the bodies' reference spheres must not overlap.

For rough surfaces, create a crater MAT file once and reuse it:

```matlab
cache_root = fullfile(pwd, 'local_cache');
if ~isfolder(cache_root), mkdir(cache_root); end
cobird.geometry.build_crater_model(cache_root, 'hemisphere_crater');
crater_file = fullfile(cache_root, 'hemisphere_crater.mat');
```

Pass this path as `crater_file` with `surface_mode="rough"`, and set
`params.roughness_f_R` to the crater-covered area fraction between 0 and 1.
The same crater model and fraction apply to both bodies.

Use `cobird.run_coupled(problem, time_step_s=120, step_count=10)` for coupled
evolution. For temperature and recoil on an existing trajectory, also supply
`state_at_time(t)` and `orbit_id` to `prepare_binary`, then call
`cobird.run_thermal` with the same time options. This driver computes temperature
and recoil along the prescribed orbit and attitude.

Both drivers return full thermal profiles in kelvin and body-frame forces [N]
and torques [N m]. Each temperature matrix has facets in rows and depth nodes
in columns, with the surface in column 1. Rough-temperature rows group crater
facets within each parent facet.
State and temperature times are in `time_s`; interval-start recoil times are
in `force_time_s`. `save_stride` controls storage, not calculation cadence.
The default keeps every step; increase it for large models.

Set `checkpoint_file="segment_01.mat"` to save a completed segment, then
`resume_file="segment_01.mat"` to continue at its endpoint. Use a new filename
for each checkpoint. Restart checks the actual crater file, body geometry,
parameters, numerical settings, and time step. `environment_id` and `orbit_id`
must identify the callback data and settings; update them if those inputs
change. No checkpoint is written unless requested.

For tides, `params.tidal.k2_over_Q` supplies a shared value; the fields
`k2_over_Q_primary` and `k2_over_Q_secondary` override it independently.
Named `RUN_TIDAL_K2_OVER_Q*` environment overrides are also supported.
The planar tidal model assumes prograde primary rotation faster than the
orbital motion.

## Interbody infrared approximation

The interbody calculation uses the center-to-center separation for dilution
and the same centerline for both projected-area cosines and self-shadow
queries, with directions transformed into the corresponding body frames.
Source emission and receiver absorption are summed separately, with cost
proportional to the total number of facets.

The approximation omits finite-body variations in pairwise distance,
direction, and visibility. Within-body radiative exchange and photon recoil
remain separate calculations. The horizon-integrated recoil formulation
follows Yan and Li (2019), *YORP effect on the evolution of asteroid spin
states*, Scientia Sinica Physica, Mechanica & Astronomica 49, 084511 (in
Chinese), [doi:10.1360/SSPMA-2019-0130](https://doi.org/10.1360/SSPMA-2019-0130).

## License and attribution

COBIRD is released under the [MIT License](LICENSE).
Copyright (c) 2026 Yixuan Wu and Jintang Li (COBIRD).
Redistributions must retain the copyright and license notices.

Third-party code and bundled PDS shape data retain their separate terms; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Software citation metadata is available in [CITATION.cff](CITATION.cff).
