# Contributing

Please open an issue before proposing a substantial change to the physical
model or output schema.

Contributions should:

- preserve SI units and state the coordinate frame of vector quantities;
- keep scientific parameters explicit and separate from numerical algorithms;
- fail explicitly when required physics or output metadata are unavailable;
- include a focused test for changed numerical behavior;
- never commit additional shape models, kernels, generated caches, compiled MEX files, or simulation results unless their redistribution and provenance have been explicitly approved; preserve the bundled PDS OBJ/XML products byte-for-byte.

## Code style

- Use `lower_snake_case` for first-party MATLAB files, functions, options, and descriptive variables.
- Scientific symbols may retain conventional capitals, with explicit body, frame, or unit suffixes such as `F_Dim_J2000` and `period_h`.
- Implement production entry points as functions and keep configuration explicit; do not clear the caller workspace or command window.
- Use MATLAB Editor smart indentation: top-level function code starts at column zero and nested blocks use four spaces. Use no tabs and keep lines within 120 characters where practical.
- Give every public function an uppercase MATLAB H1 summary and use `local_` for private helper functions.
- Use namespaced error identifiers in the form `cobird:package:function:reason` and keep comments and diagnostics in English.
- Treat persisted MAT-file fields and metadata schemas as public interfaces; validate them strictly and reject incompatible restart data.
- Preserve upstream formatting and attribution in separately licensed third-party files.

Before submitting a change, build the MEX solvers and run the full suite using
the commands in [Tests](tests/README.md).
