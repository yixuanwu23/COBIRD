%RUN_EXAMPLE Execute the official-shape pre-impact Didymos quick example.
example_root = fileparts(mfilename('fullpath'));
addpath(example_root);
result = run_didymos_pre_impact;
