function state = advance_binary_state(time_s, state, radiation, dt, model)
%ADVANCE_BINARY_STATE Integrate dynamics with interval-start radiation held fixed.
params = radiation.params;
params.et0 = time_s;
forces = model.options.forces;
if model.options.geometry == "planar"
    rhs = @(t, y) cobird.dynamics.evaluate_planar_coupled_derivative( ...
        t, y, params, forces, radiation.planar_force, radiation.moment_secondary_body(3));
else
    environment = radiation.environment;
    names = {'EMB', 'Venus', 'Jupiter'};
    positions = cell(1, 3);
    for index = 1:3
        if forces.(names{index})
            positions{index} = environment.([names{index}, '_position_m']);
        end
    end
    rhs = @(t, y) cobird.dynamics.evaluate_3d_coupled_derivative( ...
        t, y, params, forces, radiation.force_secondary_body, radiation.force_primary_inertial, ...
        radiation.moment_secondary_body, environment.primary_to_inertial(:, 3), ...
        environment.sun_position_m, positions{1}, positions{2}, positions{3});
end
if model.options.orbit_solver == "ode15s"
    [~, history] = ode15s(rhs, [0, dt], state, odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
    state = history(end, :)';
else
    count = ceil(dt/model.options.orbit_substep_s);
    h = dt/count;
    for index = 1:count
        t = (index-1)*h;
        k1 = rhs(t, state);
        k2 = rhs(t+h/2, state+h*k1/2);
        k3 = rhs(t+h/2, state+h*k2/2);
        k4 = rhs(t+h, state+h*k3);
        state = state+h*(k1+2*k2+2*k3+k4)/6;
        if model.options.geometry == "3d"
            state(7:10) = state(7:10)/norm(state(7:10));
        end
    end
end
assert(all(isfinite(state)), 'cobird:dynamics:advance_binary_state:nonfiniteState', ...
    'Dynamics produced a nonfinite state.');
end
