function result = integrate_problem(problem, mode, options)
%INTEGRATE_PROBLEM Share time stepping, output sampling, and restart I/O.
coupled = strcmp(mode, 'coupled');
required = {'initial_temperature', 'advance_temperature'};
if coupled
    required = [required, {'initial_state', 'evaluate_radiation', 'advance_state'}];
else
    required = [required, {'evaluate_forcing'}];
end
assert(all(isfield(problem, required)), 'cobird:runtime:integrate_problem:problem', ...
    'The problem definition is incomplete.');
if strlength(options.checkpoint_file) > 0
    assert(~isfile(options.checkpoint_file), 'cobird:runtime:integrate_problem:checkpointExists', ...
        'Use a new checkpoint path; existing snapshots are not overwritten.');
end
contract = [];
if isfield(problem, 'resume_contract')
    contract = problem.resume_contract;
    contract.program_id = mode;
    contract.run_config.time_step_s = options.time_step_s;
end
temperature = problem.initial_temperature;
state = [];
if coupled
    state = problem.initial_state;
end
start_time = options.start_time_s;
if strlength(options.resume_file) > 0 || strlength(options.checkpoint_file) > 0
    assert(~isempty(contract), 'cobird:runtime:integrate_problem:resumeContract', ...
        'Restart I/O requires a physical problem prepared by cobird.prepare_binary.');
end
if strlength(options.resume_file) > 0
    saved = load(options.resume_file, 'snapshot');
    snapshot = saved.snapshot;
    cobird.runtime.validate_resume_contract(snapshot.contract, contract);
    validateattributes(snapshot.time_s, {'double'}, {'scalar', 'finite', 'real'});
    start_time = snapshot.time_s;
    temperature = snapshot.temperature;
    if coupled
        validateattributes(snapshot.state, {'double'}, ...
            {'real', 'finite', 'size', size(problem.initial_state)});
        state = snapshot.state;
    end
end
if isfield(problem, 'validate_temperature')
    problem.validate_temperature(temperature);
end
save_steps = unique([0:options.save_stride:options.step_count, options.step_count]);
saved_temperature = cell(numel(save_steps), 1);
saved_state = cell(numel(save_steps), 1);
saved_temperature{1} = temperature;
saved_state{1} = state;
next_save = 2;
diagnostics = struct();
force_time = zeros(ceil(options.step_count/options.save_stride), 1);
force_primary = zeros(numel(force_time), 3);
force_secondary = force_primary;
moment_primary = force_primary;
moment_secondary = force_primary;
force_count = 0;
for step = 1:options.step_count
    time_s = start_time+(step-1)*options.time_step_s;
    if coupled
        radiation = problem.evaluate_radiation(time_s, state, temperature);
        next_state = problem.advance_state(time_s, state, radiation, options.time_step_s);
    else
        radiation = problem.evaluate_forcing(time_s, temperature);
        next_state = [];
    end
    temperature = problem.advance_temperature(time_s, temperature, radiation, options.time_step_s);
    state = next_state;
    if isfield(problem, 'summarize_radiation')
        summary = problem.summarize_radiation(radiation);
        names = fieldnames(summary);
        for index = 1:numel(names)
            name = names{index};
            if ~isfield(diagnostics, name)
                diagnostics.(name) = summary.(name);
            end
            diagnostics.(name) = max(diagnostics.(name), summary.(name));
        end
    end
    if isstruct(radiation) && isfield(radiation, 'force_primary_body') && ...
            mod(step-1, options.save_stride) == 0
        force_count = force_count+1;
        force_time(force_count) = time_s;
        force_primary(force_count, :) = radiation.force_primary_body';
        force_secondary(force_count, :) = radiation.force_secondary_body';
        moment_primary(force_count, :) = radiation.moment_primary_body';
        moment_secondary(force_count, :) = radiation.moment_secondary_body';
    end
    if next_save <= numel(save_steps) && step == save_steps(next_save)
        saved_temperature{next_save} = temperature;
        saved_state{next_save} = state;
        next_save = next_save+1;
    end
end
result = struct('time_s', start_time+save_steps(:)*options.time_step_s, ...
    'temperature', {saved_temperature}, 'final_temperature', temperature, ...
    'time_step_s', options.time_step_s, 'step_count', options.step_count, ...
    'save_stride', options.save_stride, 'diagnostics', diagnostics);
if coupled
    result.state = saved_state;
    result.final_state = state;
end
result.force_time_s = force_time(1:force_count);
result.force_primary_body_N = force_primary(1:force_count, :);
result.force_secondary_body_N = force_secondary(1:force_count, :);
result.moment_primary_body_Nm = moment_primary(1:force_count, :);
result.moment_secondary_body_Nm = moment_secondary(1:force_count, :);
snapshot = struct('time_s', start_time+options.step_count*options.time_step_s, ...
    'state', state, 'temperature', temperature, 'contract', contract);
if strlength(options.checkpoint_file) > 0
    % Save only at a completed interval boundary, never a half-updated state.
    save(options.checkpoint_file, 'snapshot', '-v7.3');
end
result.checkpoint_file = options.checkpoint_file;
end
