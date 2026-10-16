function tests = test_terrain_horizon
%TEST_TERRAIN_HORIZON Check integrated recoil against analytic horizons.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
end

function testWallBothTriangulations(test_case)
expected = local_wall_integral(2, 0.8, 1.2);
for diagonal = 1:2
    [vertices, faces] = local_wall(2, 0.8, 1.2, diagonal);
    verifyEqual(test_case, local_integrate(vertices, faces), expected, 'AbsTol', 2e-8);
end
end

function testCyclicSourceFrameAndBodyRotation(test_case)
source = [-1, -1, 0; 2, -1, 0; -1, 2, 0];
[vertices, faces] = local_wall(2, 0.8, 1.2, 1);
rotation = local_rotation([1, -2, 3], 0.83);
translation = [4, -3, 2];
source = source * rotation.' + translation;
vertices = vertices * rotation.' + translation;
expected = local_wall_integral(2, 0.8, 1.2) * rotation.';
for first = 1:3
    order = mod((first:first+2) - 1, 3) + 1;
    actual = local_body_integral(source(order, :), vertices, faces);
    verifyEqual(test_case, actual, expected, 'AbsTol', 2e-8);
end
end

function testDisjointWalls(test_case)
[first, first_faces] = local_wall(2, 0.4, 0.8, 1);
[second, second_faces] = local_wall(3, 0.6, 1.5, 2);
rotation1 = local_rotation([0, 0, 1], 0.4);
rotation2 = local_rotation([0, 0, 1], 2.3);
vertices = [first * rotation1.'; second * rotation2.'];
faces = [first_faces; second_faces + size(first, 1)];
unblocked = [0, 0, 2*pi];
expected = unblocked + ...
    (local_wall_integral(2, 0.4, 0.8) - unblocked) * rotation1.' + ...
    (local_wall_integral(3, 0.6, 1.5) - unblocked) * rotation2.';
verifyEqual(test_case, local_integrate(vertices, faces), expected, 'AbsTol', 2e-8);
end

function testEmptyAndBelowHorizon(test_case)
empty = zeros(0, 3);
verifyEqual(test_case, cobird.geometry.integrate_terrain_horizon(empty, empty, empty), ...
    [0, 0, 2*pi]);
[vertices, faces] = local_wall(2, 0.8, 1.2, 1);
vertices(:, 3) = vertices(:, 3) - 1.2;
verifyEqual(test_case, local_integrate(vertices, faces), [0, 0, 2*pi], 'AbsTol', 1e-12);
end

function testNarrowWall(test_case)
[vertices, faces] = local_wall(1, 1e-5, 2, 1);
rotation = local_rotation([0, 0, 1], 0.731);
expected = local_wall_integral(1, 1e-5, 2) * rotation.';
verifyEqual(test_case, local_integrate(vertices * rotation.', faces), ...
    expected, 'AbsTol', 2e-10);
end

function testClippedTriangle(test_case)
D = 2; W = 0.8; H = 1.2;
vertices = [D, -W, -H; D, W, H; D, W, -H];
% Above z=0, this wall has h(phi)=atan((H/W)*sin(phi)), 0<phi<atan(W/D).
alpha = atan(W/D);
q = H/W;
s = sin(alpha); c = cos(alpha);
A = 1 + q^2;
r = sqrt(1 + q^2*s^2);
expected = [s/r - s, 1/A - c/(A*r) - (1-c), ...
    2*pi - asin(q/sqrt(A)) + asin(q*c/sqrt(A)) + q/A - q*c/(A*r)];
rotation = local_rotation([0, 0, 1], 0.37);
verifyEqual(test_case, local_integrate(vertices * rotation.', [1, 2, 3]), ...
    expected * rotation.', 'AbsTol', 2e-8);
end

function testOverheadTriangle(test_case)
% The upward axis crosses the triangle interior, so every azimuth reaches pi/2.
vertices = [-2, -1, 1; 2, -1, 1; 0, 2, 1];
verifyEqual(test_case, local_integrate(vertices, [1, 2, 3]), [0, 0, 0], 'AbsTol', 2e-8);
end

function testVerticalPlaneHasZeroAngularMeasure(test_case)
% This triangle crosses the upward axis only in the x-z plane.
vertices = [-1, 0, 1; 1, 0, 1; 0, 0, 2];
verifyEqual(test_case, local_integrate(vertices, [1, 2, 3]), ...
    [0, 0, 2*pi], 'AbsTol', 2e-8);
vertices = [1, 2, 1; -3, -6, 1; 2, 4, 2];
verifyEqual(test_case, local_integrate(vertices, [1, 2, 3]), ...
    [0, 0, 2*pi], 'AbsTol', 2e-8);
end

function testUpwardAxisOnTriangleEdge(test_case)
vertices = [-1, 0, 1; 1, 0, 1; 0, 2, 1];
verifyEqual(test_case, local_integrate(vertices, [1, 2, 3]), ...
    [0, -2, pi], 'AbsTol', 2e-8);
end

function [vertices, faces] = local_wall(D, W, H, diagonal)
vertices = [D, -W, 0; D, W, 0; D, W, H; D, -W, H];
if diagonal == 1
    faces = [1, 2, 3; 1, 3, 4];
else
    faces = [1, 2, 4; 2, 3, 4];
end
end

function iv = local_wall_integral(D, W, H)
alpha = atan(W/D);
k = H/D;
s = sin(alpha);
C = (1+k^2) * sqrt(1+k^2*cos(alpha)^2);
iv = [2*s*(1/C-1), 0, 2*pi - 2*asin(k*s/sqrt(1+k^2)) + 2*k*s/C];
end

function rotation = local_rotation(axis, angle)
axis = axis / norm(axis);
skew = [0, -axis(3), axis(2); axis(3), 0, -axis(1); -axis(2), axis(1), 0];
rotation = cos(angle)*eye(3) + (1-cos(angle))*(axis.'*axis) + sin(angle)*skew;
end

function iv = local_integrate(vertices, faces)
iv = cobird.geometry.integrate_terrain_horizon(vertices(faces(:, 1), :), ...
    vertices(faces(:, 2), :), vertices(faces(:, 3), :));
end

function iv = local_body_integral(source, vertices, faces)
edge = source(2, :) - source(1, :);
normal = cross(edge, source(3, :) - source(1, :));
normal = normal / norm(normal);
x_local = edge / norm(edge);
frame = [x_local; cross(normal, x_local); normal].';
local_vertices = (vertices - mean(source, 1)) * frame;
iv = local_integrate(local_vertices, faces) * frame.';
end
