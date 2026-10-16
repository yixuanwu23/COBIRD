function tests = test_mutual_visibility
%TEST_MUTUAL_VISIBILITY Test mutual-radiation geometry.
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'src'));
test_case.TestData.solve = @cobird.geometry.compute_mutual_visibility;
end

function testAnalyticAndSmallPositiveCosines(test_case)
primary_normals = [ ...
    1, 0, 0; ...
    0.01, sqrt(1 - 0.01^2), 0; ...
    0, 1, 0; ...
   -1, 0, 0];
secondary_normals = [-1, 0, 0];

[primary_indices, secondary_indices, primary_cosines, ...
    secondary_cosines, distance] = solve_visibility( ...
    test_case, primary_normals, secondary_normals, [10; 0; 0]);

verifyEqual(test_case, primary_indices, uint32([1; 2]));
verifyEqual(test_case, secondary_indices, uint32(1));
verifyEqual(test_case, primary_cosines, [1; 0.01], 'AbsTol', 1e-15);
verifyEqual(test_case, secondary_cosines, 1);
verifyEqual(test_case, distance, 10);

projected_sum = secondary_cosines * ...
    sum(primary_cosines .* [3; 7]) / (pi * distance^2);
verifyEqual(test_case, projected_sum, 3.07 / (100 * pi), 'AbsTol', 1e-15);
end

function testCentroidAndRadiusIndependence(test_case)
primary_normals = [0.6, 0.8, 0];
secondary_normals = -primary_normals;
separation = [10; 1; 2];
[primary_reabsorb, primary_shadow] = create_clear_shadow(1);
[secondary_reabsorb, secondary_shadow] = create_clear_shadow(1);
solve = test_case.TestData.solve;

baseline = cell(1, 5);
[baseline{:}] = solve( ...
    [0, 0, 0], primary_normals, primary_reabsorb, ...
    [10, 1, 2], secondary_normals, secondary_reabsorb, ...
    separation, 1, 2, primary_shadow, secondary_shadow, ...
    -separation / norm(separation));

changed = cell(1, 5);
[changed{:}] = solve( ...
    [8, -100, 5], primary_normals, primary_reabsorb, ...
    [-20, 300, 2], secondary_normals, secondary_reabsorb, ...
    separation, 1000, 2000, primary_shadow, secondary_shadow, ...
    -separation / norm(separation));

verifyEqual(test_case, changed, baseline);
end

function testInverseSquareAndReciprocity(test_case)
[~, ~, primary_cosine, secondary_cosine, distance] = solve_visibility( ...
    test_case, [1, 0, 0], [-1, 0, 0], [10; 0; 0]);
[~, ~, primary_cosine_2, secondary_cosine_2, distance_2] = ...
    solve_visibility(test_case, [1, 0, 0], [-1, 0, 0], [20; 0; 0]);

exchange = primary_cosine * secondary_cosine / distance^2;
exchange_2 = primary_cosine_2 * secondary_cosine_2 / distance_2^2;
verifyEqual(test_case, exchange_2, exchange / 4, 'AbsTol', 1e-15);

[~, ~, reverse_primary_cosine, reverse_secondary_cosine, reverse_distance] = ...
    solve_visibility(test_case, [-1, 0, 0], [1, 0, 0], [-10; 0; 0]);
primary_area = 2;
secondary_area = 7;
forward = primary_area * secondary_area * primary_cosine * secondary_cosine / ...
    (pi * distance^2);
reverse = secondary_area * primary_area * reverse_primary_cosine * ...
    reverse_secondary_cosine / (pi * reverse_distance^2);
verifyEqual(test_case, reverse, forward, 'AbsTol', 1e-15);
end

function testSeparateBodyFrameShadowDirections(test_case)
direction = [0.8; 0.5; 0.2];
direction = direction / norm(direction);
rotation_angle = 0.71;
secondary_to_primary = [ ...
    cos(rotation_angle), -sin(rotation_angle), 0; ...
    sin(rotation_angle),  cos(rotation_angle), 0; ...
    0,                    0,                   1];
primary_normals = repmat(direction', 3, 1);
secondary_normals_in_primary = -primary_normals;

[primary_reabsorb, primary_shadow] = create_clear_shadow(3);
[secondary_reabsorb, secondary_shadow] = create_clear_shadow(3);
primary_reabsorb.has_candidates = true(3, 1);
secondary_reabsorb.has_candidates = true(3, 1);
primary_shadow.occ_rev_map = uint32((1:3)');
secondary_shadow.occ_rev_map = uint32((1:3)');
primary_shadow.shadow_table = true(3, 32);
secondary_shadow.shadow_table = true(3, 32);
primary_shadow.shadow_table([1, 3], ...
    direction_bin(primary_shadow, direction)) = false;
direction_to_primary_secondary = secondary_to_primary' * (-direction);
secondary_shadow.shadow_table([2, 3], ...
    direction_bin(secondary_shadow, direction_to_primary_secondary)) = false;

[primary_indices, secondary_indices, primary_cosines, secondary_cosines] = ...
    test_case.TestData.solve( ...
    zeros(3), primary_normals, primary_reabsorb, ...
    zeros(3), secondary_normals_in_primary, secondary_reabsorb, ...
    10 * direction, 1, 1, primary_shadow, secondary_shadow, ...
    direction_to_primary_secondary);

verifyEqual(test_case, primary_indices, uint32([1; 3]));
verifyEqual(test_case, secondary_indices, uint32([2; 3]));
verifyEqual(test_case, primary_cosines, ones(2, 1), 'AbsTol', 1e-14);
verifyEqual(test_case, secondary_cosines, ones(2, 1), 'AbsTol', 1e-14);
end

function testRigidFrameRotation(test_case)
primary_normals = [0.8, 0.6, 0; 0.4, 0.3, sqrt(0.75)];
secondary_normals = -primary_normals;
separation = [10; 2; 1];

baseline = cell(1, 5);
[baseline{:}] = solve_visibility( ...
    test_case, primary_normals, secondary_normals, separation);

axis = [0.2; 0.3; 0.9];
axis = axis / norm(axis);
cross_matrix = [ ...
    0,       -axis(3),  axis(2); ...
    axis(3),  0,       -axis(1); ...
   -axis(2),  axis(1),  0];
rotation = eye(3) + sin(0.6) * cross_matrix + ...
    (1 - cos(0.6)) * cross_matrix^2;

rotated = cell(1, 5);
[rotated{:}] = solve_visibility( ...
    test_case, primary_normals * rotation', secondary_normals * rotation', ...
    rotation * separation);

verifyEqual(test_case, rotated{1}, baseline{1});
verifyEqual(test_case, rotated{2}, baseline{2});
for output_index = 3:5
    verifyEqual(test_case, rotated{output_index}, baseline{output_index}, ...
        'AbsTol', 2e-14);
end
end

function testEmptyFacingSet(test_case)
outputs = cell(1, 5);
[outputs{:}] = solve_visibility( ...
    test_case, [-1, 0, 0], [-1, 0, 0], [10; 0; 0]);
verifyTrue(test_case, all(cellfun(@isempty, outputs)));
end

function testEmptyShadowedSet(test_case)
[primary_reabsorb, primary_shadow] = create_clear_shadow(1);
[secondary_reabsorb, secondary_shadow] = create_clear_shadow(1);
primary_reabsorb.has_candidates = true;
primary_shadow.occ_rev_map = uint32(1);
primary_shadow.shadow_table = true(1, 32);

outputs = cell(1, 5);
[outputs{:}] = test_case.TestData.solve( ...
    zeros(1, 3), [1, 0, 0], primary_reabsorb, ...
    zeros(1, 3), [-1, 0, 0], secondary_reabsorb, ...
    [10; 0; 0], 1, 1, primary_shadow, secondary_shadow, [-1; 0; 0]);
verifyTrue(test_case, all(cellfun(@isempty, outputs)));
end

function testInvalidSeparation(test_case)
verifyError(test_case, ...
    @() solve_visibility(test_case, [1, 0, 0], [-1, 0, 0], [0; 0; 0]), ...
    'cobird:geometry:compute_mutual_visibility:invalidSeparation');
verifyError(test_case, ...
    @() solve_visibility(test_case, [1, 0, 0], [-1, 0, 0], [NaN; 0; 0]), ...
    'cobird:geometry:compute_mutual_visibility:invalidSeparation');
end

function testInvalidShadowMapping(test_case)
[primary_reabsorb, primary_shadow] = create_clear_shadow(1);
[secondary_reabsorb, secondary_shadow] = create_clear_shadow(1);
primary_reabsorb.has_candidates = true;

verifyError(test_case, @() test_case.TestData.solve( ...
    zeros(1, 3), [1, 0, 0], primary_reabsorb, ...
    zeros(1, 3), [-1, 0, 0], secondary_reabsorb, ...
    [10; 0; 0], 1, 1, primary_shadow, secondary_shadow, [-1; 0; 0]), ...
    'cobird:geometry:compute_mutual_visibility:invalidShadowMap');
end

function testInvalidSecondaryDirection(test_case)
[primary_reabsorb, primary_shadow] = create_clear_shadow(1);
[secondary_reabsorb, secondary_shadow] = create_clear_shadow(1);

verifyError(test_case, @() test_case.TestData.solve( ...
    zeros(1, 3), [1, 0, 0], primary_reabsorb, ...
    zeros(1, 3), [-1, 0, 0], secondary_reabsorb, ...
    [10; 0; 0], 1, 1, primary_shadow, secondary_shadow, [0; 0; 0]), ...
    'cobird:geometry:compute_mutual_visibility:invalidDimorphosDirection');
end

function varargout = solve_visibility( ...
    test_case, primary_normals, secondary_normals, separation)
[primary_reabsorb, primary_shadow] = create_clear_shadow(size(primary_normals, 1));
[secondary_reabsorb, secondary_shadow] = ...
    create_clear_shadow(size(secondary_normals, 1));
direction_to_primary_secondary = -separation / norm(separation);

if nargout == 0
    test_case.TestData.solve( ...
        zeros(size(primary_normals)), primary_normals, primary_reabsorb, ...
        zeros(size(secondary_normals)), secondary_normals, secondary_reabsorb, ...
        separation, 1, 1, primary_shadow, secondary_shadow, ...
        direction_to_primary_secondary);
else
    [varargout{1:nargout}] = test_case.TestData.solve( ...
        zeros(size(primary_normals)), primary_normals, primary_reabsorb, ...
        zeros(size(secondary_normals)), secondary_normals, secondary_reabsorb, ...
        separation, 1, 1, primary_shadow, secondary_shadow, ...
        direction_to_primary_secondary);
end
end

function [reabsorb, shadow] = create_clear_shadow(facet_count)
reabsorb = struct('has_candidates', false(facet_count, 1));
shadow = struct( ...
    'N_az', 8, ...
    'N_el', 4, ...
    'daz', 2 * pi / 8, ...
    'del', pi / 4, ...
    'occ_rev_map', zeros(facet_count, 1, 'uint32'), ...
    'shadow_table', false(facet_count, 32));
end

function bin = direction_bin(shadow, direction)
azimuth = atan2(direction(2), direction(1));
elevation = asin(max(-1, min(1, direction(3))));
azimuth_bin = max(1, min(shadow.N_az, ...
    floor((azimuth + pi) / shadow.daz) + 1));
elevation_bin = max(1, min(shadow.N_el, ...
    floor((elevation + pi/2) / shadow.del) + 1));
bin = (elevation_bin - 1) * shadow.N_az + azimuth_bin;
end
