function tests = test_minimal_binary
%TEST_MINIMAL_BINARY Exercise the synthetic continuous-coupling workflow.
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
test_root = fileparts(mfilename('fullpath'));
project_root = fileparts(test_root);
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'examples', 'minimal_binary'));
test_case.TestData.project_root = project_root;
end

function testSyntheticEndToEndRun(test_case)
result = run_minimal_binary(step_count=24, verbose=false);
check_values = struct2cell(result.checks);
verifyTrue(test_case, all([check_values{:}]));
verifyEqual(test_case, result.schema, 'COBIRD_MINIMAL_BINARY');
verifyEqual(test_case, result.mutual_geometry, 'common_centerline');
verifyGreaterThan(test_case, result.maximum_shadowed_faces, 0);
verifyGreaterThan(test_case, result.maximum_mutual_flux_W_m2, 0);
verifyGreaterThan(test_case, result.maximum_rough_recoil_difference_N, 0);
verifyGreaterThan(test_case, result.surface_temperature_range_K, 1);
end
