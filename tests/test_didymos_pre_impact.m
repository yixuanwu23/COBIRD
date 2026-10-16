function tests = test_didymos_pre_impact
%TEST_DIDYMOS_PRE_IMPACT Verify PDS data integrity and the public example.
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
test_root = fileparts(mfilename('fullpath'));
project_root = fileparts(test_root);
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'examples', 'didymos_pre_impact'));
test_case.TestData.project_root = project_root;
end

function testPdsSourceIntegrity(test_case)
data_root = fullfile(test_case.TestData.project_root, ...
    'data', 'shape_models', 'dart_pre_impact');
records = { ...
    'didymos_g_9309mm_spc_obj_0000n00000_v003.obj', ...
    'b5349374186d254ca21b7b00ddd9f6dea7e2762ca6bbff13f43151ee448e3f20'; ...
    'didymos_g_9309mm_spc_obj_0000n00000_v003.xml', ...
    'f09a8a3364c6fc9e9c6df54206b03e550f0da16a896d3d65c6e5d74b0d525564'; ...
    'dimorphos_g_1940mm_spc_obj_0000n00000_v004.obj', ...
    '833eaeac8c49338abe64012a679b1b9695f774f1b5e4842a8e5b9f0af1d28664'; ...
    'dimorphos_g_1940mm_spc_obj_0000n00000_v004.xml', ...
    '715a4f105d2284ad818b5a3b5eb64087be8f37a6c997fafc31e88724a2c765f6'};
for record_index = 1:size(records, 1)
    path = fullfile(data_root, records{record_index, 1});
    verifyEqual(test_case, cobird.io.compute_file_sha256(path), ...
        records{record_index, 2});
end

end

function testOfficialShapeExample(test_case)
result = run_didymos_pre_impact( ...
    step_count=2, time_step_s=300, verbose=false);
check_values = struct2cell(result.checks);
verifyTrue(test_case, all([check_values{:}]));
verifyEqual(test_case, result.schema, 'COBIRD_DIDYMOS_PRE_IMPACT');
verifyEqual(test_case, result.mutual_geometry, 'common_centerline');
verifyEqual(test_case, result.didymos.vertex_count, 24578);
verifyEqual(test_case, result.didymos.face_count, 49152);
verifyEqual(test_case, result.dimorphos.vertex_count, 24578);
verifyEqual(test_case, result.dimorphos.face_count, 49152);
end
