# NAIF MICE

MICE is an optional external dependency for SPICE-driven workflows and is not
distributed with COBIRD. Download the platform-specific MATLAB Toolkit from
the [official NAIF page](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html).

Point COBIRD to the extracted `mice` directory using either:

```matlab
startup_cobird(mice_root="/absolute/path/to/mice")
```

or the `COBIRD_MICE_ROOT` environment variable. The selected root must contain
`src/mice/` and `lib/mice.<mexext>`.

For a repository-local installation, use one of these ignored directories:

```text
external/mice/win64/mice/
external/mice/linux64/mice/
external/mice/mac64/mice/
```

COBIRD never downloads or executes MICE automatically.
