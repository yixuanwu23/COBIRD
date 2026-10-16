function result = run_coupled(problem, options)
%RUN_COUPLED Advance synchronized thermal and orbit-attitude states.
% Build a physical problem with cobird.prepare_binary. Radiation uses the
% interval-start state; saved temperatures and states share endpoint times.
% checkpoint_file saves the completed run; resume_file continues from it.

arguments
    problem (1,1) struct
    options.time_step_s (1,1) double {mustBeFinite,mustBePositive}
    options.step_count (1,1) double {mustBeInteger,mustBeNonnegative}
    options.save_stride (1,1) double {mustBeInteger,mustBePositive} = 1
    options.start_time_s (1,1) double {mustBeFinite} = 0
    options.resume_file (1,1) string = ""
    options.checkpoint_file (1,1) string = ""
end
result = cobird.runtime.integrate_problem(problem,'coupled',options);
end
