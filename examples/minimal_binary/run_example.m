%RUN_EXAMPLE Execute the data-free COBIRD quick-start example.
example_root = fileparts(mfilename('fullpath'));
addpath(example_root);
result = run_minimal_binary;
