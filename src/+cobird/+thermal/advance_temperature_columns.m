function [result, surface_flux] = advance_temperature_columns( ...
        ~, mu_0, T_all_matrix, params, d_sun_AU, dt, Q_terrain, ...
        dz_ov, Nz_ov, solver_config)
%ADVANCE_TEMPERATURE_COLUMNS Advance independent one-dimensional thermal columns.
% Implicit mode uses backward Euler with a linearized Stefan-Boltzmann
% surface boundary. Without a time step, explicit mode returns dT/dt.
% surface_flux is thermal emission [W/m^2] at the updated temperature in
% implicit mode, or at the input temperature in explicit mode.
% Q_terrain is additional absorbed flux [W/m^2]. dz_ov may be a uniform
% spacing or a depth-node vector beginning at z=0.

persistent has_fixed_solver has_generic_solver has_nonuniform_solver
persistent grid_cache_z_bank grid_cache_interval_bank

N_faces = size(T_all_matrix, 1);
implicit_mode = (nargin >= 6) && ~isempty(dt);

if nargin < 7 || isempty(Q_terrain)
    Q_terrain = zeros(N_faces, 1);
end
if nargin < 10 || isempty(solver_config)
    solver_config = struct('thread_count', 1, 'tile_size', 16);
end

if isfield(params, 'use_roughness_A_correction')
    error('cobird:thermal:advance_temperature_columns:unsupportedRoughnessCorrection', ...
        ['use_roughness_A_correction is not supported. Select surface_mode=''rough'' ' ...
         'for explicit surface roughness.']);
end

q = (1 - params.A) * params.F_s / d_sun_AU^2;

k_therm = params.k;
rho = params.rho;
Cp = params.Cp;
if nargin >= 9 && ~isempty(Nz_ov)
    Nz = Nz_ov;
else
    Nz = params.Nz;
end
if nargin >= 8 && ~isempty(dz_ov)
    dz = dz_ov;
else
    dz = params.dz;
end

nonuniform_grid = ~isscalar(dz);
if nonuniform_grid
    z = dz(:);
    if numel(z) ~= Nz || Nz < 2 || z(1) ~= 0 || ...
            any(~isfinite(z)) || any(diff(z) <= 0)
        error('cobird:thermal:advance_temperature_columns:badDepthGrid', ...
            ['A nonuniform depth grid must be a finite, strictly increasing ' ...
             'Nz-element vector that begins at zero.']);
    end
end

eps_sig = params.epsilon * params.sigma;
E_in = q .* mu_0 + Q_terrain;

if implicit_mode
    if nonuniform_grid
        alpha = k_therm / (rho * Cp);
        rho_cp = rho * Cp;
        if isempty(grid_cache_z_bank)
            grid_cache_z_bank = {};
            grid_cache_interval_bank = {};
        end
        grid_slot = [];
        for candidate = 1:numel(grid_cache_z_bank)
            if isequal(grid_cache_z_bank{candidate}, z)
                grid_slot = candidate;
                break
            end
        end
        if isempty(grid_slot)
            interval = diff(z);
            if numel(grid_cache_z_bank) == 2
                grid_cache_z_bank(1) = [];
                grid_cache_interval_bank(1) = [];
            end
            grid_slot = numel(grid_cache_z_bank) + 1;
            grid_cache_z_bank{grid_slot} = z;
            grid_cache_interval_bank{grid_slot} = interval;
        else
            interval = grid_cache_interval_bank{grid_slot};
        end
        [lower_coeff, upper_coeff, surface_beta] = ...
            local_nonuniform_coefficients(interval, alpha, rho_cp);

        if isempty(has_nonuniform_solver) || ~has_nonuniform_solver
            has_nonuniform_solver = ~isempty(which( ...
                'cobird.thermal.calculate_temperature_nonuniform_threads_mex'));
        end
        inv_dt = 1.0 / dt;
        if ~has_nonuniform_solver
            error('cobird:thermal:advance_temperature_columns:missingMex', ...
                ['The nonuniform-grid thermal-column MEX solver is missing. ' ...
                 'Run cobird.thermal.build_mex_solvers first.']);
        end
        if nargout > 1
            [result, surface_flux] = ...
                cobird.thermal.calculate_temperature_nonuniform_threads_mex( ...
                T_all_matrix, mu_0, Q_terrain, q, eps_sig, ...
                lower_coeff, upper_coeff, surface_beta, inv_dt, ...
                solver_config.tile_size, solver_config.thread_count);
        else
            result = cobird.thermal.calculate_temperature_nonuniform_threads_mex( ...
                T_all_matrix, mu_0, Q_terrain, q, eps_sig, ...
                lower_coeff, upper_coeff, surface_beta, inv_dt, ...
                solver_config.tile_size, solver_config.thread_count);
        end
        return
    end

    kappa = k_therm / (rho * Cp * dz^2);
    beta = 1.0 / (rho * Cp * dz);
    inv_dt = 1.0 / dt;

    if isempty(has_fixed_solver) || ~has_fixed_solver
        has_fixed_solver = ~isempty(which('cobird.thermal.calculate_temperature_fixed20_threads_mex'));
    end
    if isempty(has_generic_solver) || ~has_generic_solver
        has_generic_solver = ~isempty(which('cobird.thermal.calculate_temperature_tiled_threads_mex'));
    end
    thread_count = solver_config.thread_count;
    tile_size = solver_config.tile_size;
    use_fixed_solver = (Nz == 20) && (tile_size == 16);
    if use_fixed_solver && has_fixed_solver
        if nargout > 1
            [result, surface_flux] = cobird.thermal.calculate_temperature_fixed20_threads_mex( ...
                T_all_matrix, mu_0, Q_terrain, q, eps_sig, kappa, beta, inv_dt, tile_size, ...
                thread_count);
        else
            result = cobird.thermal.calculate_temperature_fixed20_threads_mex(T_all_matrix, mu_0, Q_terrain, ...
                q, eps_sig, kappa, beta, inv_dt, tile_size, ...
                thread_count);
        end
        return
    end
    if ~has_generic_solver
        error('cobird:thermal:advance_temperature_columns:missingMex', ...
            ['The general thermal-column MEX solver is missing. ' ...
             'Run cobird.thermal.build_mex_solvers first.']);
    end
    if nargout > 1
        [result, surface_flux] = cobird.thermal.calculate_temperature_tiled_threads_mex( ...
            T_all_matrix, mu_0, Q_terrain, q, eps_sig, kappa, beta, inv_dt, ...
            tile_size, thread_count);
    else
        result = cobird.thermal.calculate_temperature_tiled_threads_mex( ...
            T_all_matrix, mu_0, Q_terrain, q, eps_sig, kappa, beta, inv_dt, ...
            tile_size, thread_count);
    end
    return
end

if nonuniform_grid
    error('cobird:thermal:advance_temperature_columns:nonuniformExplicitUnsupported', ...
        'Nonuniform depth grids support only implicit time integration.');
end

dT_dt = zeros(N_faces, Nz);
if Nz == 1
    dT_dt(:,1) = (E_in - eps_sig .* T_all_matrix(:,1).^4) / (rho * Cp * dz);
else
    coeff = k_therm / (rho * Cp * dz^2);
    dT_dt(:,1) = (k_therm / dz * (T_all_matrix(:,2) - T_all_matrix(:,1)) + ...
                  E_in - eps_sig .* T_all_matrix(:,1).^4) / (rho * Cp * dz);
    dT_dt(:,2:Nz-1) = coeff * (T_all_matrix(:,3:Nz) ...
        - 2*T_all_matrix(:,2:Nz-1) + T_all_matrix(:,1:Nz-2));
    dT_dt(:,Nz) = coeff * (T_all_matrix(:,Nz-1) - T_all_matrix(:,Nz));
end
result = dT_dt;
if nargout > 1
    surface_flux = eps_sig .* T_all_matrix(:,1).^4;
end
end

function [lower_coeff, upper_coeff, surface_beta] = ...
        local_nonuniform_coefficients(interval, alpha, rho_cp)
Nz = numel(interval) + 1;
lower_coeff = zeros(Nz, 1);
upper_coeff = zeros(Nz, 1);

% Match the uniform-grid surface control layer and insulated bottom.
upper_coeff(1) = alpha / interval(1)^2;
surface_beta = 1 / (rho_cp * interval(1));

if Nz > 2
    left = interval(1:end-1);
    right = interval(2:end);
    lower_coeff(2:Nz-1) = 2 * alpha ./ (left .* (left + right));
    upper_coeff(2:Nz-1) = 2 * alpha ./ (right .* (left + right));
end
lower_coeff(Nz) = alpha / interval(end)^2;
end
