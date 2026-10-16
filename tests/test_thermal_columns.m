function tests = test_thermal_columns
%TEST_THERMAL_COLUMNS Check column discretization and compiled input buffers.
tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
end

function testChangingMaterialAndDepthGrids(test_case)
params = struct('A', 0.1, 'F_s', 1361, 'epsilon', 0.9, ...
    'sigma', 5.670374419e-8, 'k', 0.05, 'rho', 1800, 'Cp', 680);
grids = {[0; 0.02; 0.05; 0.1], [0; 0.03; 0.08; 0.15], ...
    [0; 0.01; 0.04; 0.12]};
sequence = [1, 2, 1, 3, 2, 1];
temperature = [250, 240, 230, 225; 260, 245, 235, 228];
mu = [0.8; 0.2];
extra = [10; 3];
distance_au = 1.2;
dt = 120;
for step = 1:numel(sequence)
    z = grids{sequence(step)};
    params.k = 0.04 + 0.01*step;
    params.rho = 1800 + 20*step;
    params.Cp = 680 + 2*step;
    [actual, emission] = cobird.thermal.advance_temperature_columns( ...
        [], mu, temperature, params, distance_au, dt, extra, z, numel(z));
    expected = local_column_reference(temperature, mu, extra, params, z, distance_au, dt);
    verifyEqual(test_case, actual, expected, 'AbsTol', 1e-10);
    verifyEqual(test_case, emission, params.epsilon*params.sigma*actual(:, 1).^4, ...
        'RelTol', 1e-14);
end
end

function testExplicitDerivativeAndEmission(test_case)
params = struct('A', 0.1, 'F_s', 1361, 'epsilon', 0.9, ...
    'sigma', 5.670374419e-8, 'k', 0.05, 'rho', 1800, 'Cp', 680, 'dz', 0.02);
mu = [0.8; 0.2];
extra = [10; 3];
distance_au = 1.2;
for depth_count = [1, 2, 4, 20]
    params.Nz = depth_count;
    temperature = [250; 260] - 2 * (0:depth_count-1);
    [derivative, emission] = cobird.thermal.advance_temperature_columns( ...
        [], mu, temperature, params, distance_au, [], extra);
    expected_emission = params.epsilon * params.sigma * temperature(:, 1).^4;
    absorbed = (1 - params.A) * params.F_s / distance_au^2 * mu + extra;
    expected = zeros(size(temperature));
    expected(:, 1) = (absorbed - expected_emission) / (params.rho * params.Cp * params.dz);
    if depth_count > 1
        coefficient = params.k / (params.rho * params.Cp * params.dz^2);
        expected(:, 1) = expected(:, 1) + coefficient * (temperature(:, 2) - temperature(:, 1));
        expected(:, 2:end-1) = coefficient * (temperature(:, 3:end) ...
            - 2 * temperature(:, 2:end-1) + temperature(:, 1:end-2));
        expected(:, end) = coefficient * (temperature(:, end-1) - temperature(:, end));
    end
    verifyEqual(test_case, derivative, expected, 'AbsTol', 1e-12);
    verifyEqual(test_case, emission, expected_emission, 'RelTol', 1e-14);
end
end

function testThermalMexRejectsInvalidBuffers(test_case)
names = {'calculate_temperature_tiled_threads_mex', ...
    'calculate_temperature_fixed20_threads_mex', ...
    'calculate_temperature_nonuniform_threads_mex'};
temperature = reshape(linspace(190, 250, 60), 3, 20);
uniform = {temperature, [1; 0.5; 0], [10; 2; 0], 1200, ...
    0.9*5.670374419e-8, 1e-5, 1e-4, 0.1, 16, 1};
nonuniform = [uniform(1:5), ...
    {[0; 1e-5*ones(19, 1)], [1e-5*ones(19, 1); 0], 1e-4, 0.1, 16, 1}];
conversions = {@sparse, @single, @(value) complex(value, ones(size(value)))};
for kernel = 1:numel(names)
    solver = str2func(['cobird.thermal.', names{kernel}]);
    prefix = ['cobird:thermal:', names{kernel}, ':'];
    inputs = uniform;
    positions = [1, 2, 3];
    identifiers = {'temperature', 'muSize', 'fluxSize'};
    if kernel == 3
        inputs = nonuniform;
        positions = [positions, 6, 7];
        identifiers = [identifiers, {'lowerSize', 'upperSize'}];
    end
    for buffer = 1:numel(positions)
        position = positions(buffer);
        for conversion = 1:numel(conversions)
            invalid = inputs;
            invalid{position} = conversions{conversion}(inputs{position});
            verifyError(test_case, @() solver(invalid{:}), [prefix, identifiers{buffer}]);
        end
    end
    invalid = inputs;
    invalid{1} = reshape(temperature, 3, 4, 5);
    verifyError(test_case, @() solver(invalid{:}), [prefix, 'temperature']);
end
end

function testGeneralThermalMexRejectsZeroDepth(test_case)
verifyError(test_case, @() cobird.thermal.calculate_temperature_tiled_threads_mex( ...
    zeros(3, 0), ones(3, 1), zeros(3, 1), 1200, 0.9*5.670374419e-8, ...
    1e-5, 1e-4, 0.1, 16, 1), ...
    'cobird:thermal:calculate_temperature_tiled_threads_mex:temperature');
end

function result = local_column_reference(temperature, mu, extra, params, z, distance_au, dt)
spacing = diff(z);
capacity = params.rho*params.Cp * ...
    [spacing(1); (spacing(1:end-1) + spacing(2:end))/2; spacing(end)];
conductance = params.k ./ spacing;
matrix = diag(capacity/dt + [conductance; 0] + [0; conductance]) ...
    - diag(conductance, 1) - diag(conductance, -1);
absorbed = (1 - params.A)*params.F_s/distance_au^2*mu + extra;
epsilon_sigma = params.epsilon*params.sigma;
result = zeros(size(temperature));
for facet = 1:size(temperature, 1)
    system = matrix;
    surface = temperature(facet, 1);
    system(1, 1) = system(1, 1) + 4*epsilon_sigma*surface^3;
    rhs = capacity/dt .* temperature(facet, :)';
    rhs(1) = rhs(1) + absorbed(facet) + 3*epsilon_sigma*surface^4;
    result(facet, :) = (system\rhs)';
end
end
