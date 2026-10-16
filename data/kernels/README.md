# SPICE kernels

COBIRD does not redistribute SPICE kernels. Download public kernels from the
[NAIF archive](https://naif.jpl.nasa.gov/naif/data.html), retain their original
filenames and metadata, and keep them outside Git.

For a SPICE-based case, copy `cobird.tm.template` to an untracked `.tm` file,
replace every placeholder, and use absolute paths. Record the source URL,
solution/version, coverage interval, and SHA-256 value in a local copy of
`kernel_manifest.example.json`.

The manifest is a manual provenance record. COBIRD does not read it or verify
kernel coverage or checksums at runtime. Select kernels that cover the full
requested interval and include the required targets and frames.

After configuring [MICE](../../external/mice/README.md), load the meta-kernel:

```matlab
cobird.spice.load_meta_kernel("/absolute/path/to/case.tm")
```

Supply SPICE results through `environment(t)` in `cobird.prepare_binary`.
Define how simulation seconds map to SPICE ephemeris time, convert kilometre
positions to metres, and request primary-relative Sun and perturber vectors.
Use one inertial frame consistently for these vectors, the orbital state,
and body-to-inertial rotations. The core does not select a default frame.
The optional `cobird.attitude.compute_didymos_spin_vector` helper returns
vectors in `ECLIPJ2000`.
Include the kernel selection, epoch, and frame in `environment_id` so restart
checks can distinguish callback configurations.

Both bundled examples and the test suite run without SPICE kernels.
