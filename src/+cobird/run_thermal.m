function result = run_thermal(problem, options)
%RUN_THERMAL Advance temperatures along a prescribed orbit and attitude.
% Build a physical problem with cobird.prepare_binary and state_at_time.
% Every interval advances temperature; save_stride controls output only.
% Recoil is reported but never changes the prescribed trajectory.

arguments
    problem (1,1) struct
    options.time_step_s (1,1) double {mustBeFinite,mustBePositive}
    options.step_count (1,1) double {mustBeInteger,mustBeNonnegative}
    options.save_stride (1,1) double {mustBeInteger,mustBePositive} = 1
    options.start_time_s (1,1) double {mustBeFinite} = 0
    options.resume_file (1,1) string = ""
    options.checkpoint_file (1,1) string = ""
end
result = cobird.runtime.integrate_problem(problem,'thermal',options);
end
