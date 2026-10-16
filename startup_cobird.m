function project_root = startup_cobird(options)
%STARTUP_COBIRD Add COBIRD to the MATLAB path and detect optional MICE.

arguments
    options.mice_root (1,1) string = ""
    options.quiet (1,1) logical = false
end

project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'setup'));

mice_root = cobird.runtime.setup_project_paths( ...
    project_root, require_mice=false, mice_root=char(options.mice_root));
if ~options.quiet
    if isempty(mice_root)
        fprintf(['COBIRD is ready. MICE is optional for the bundled examples ' ...
            'and is not currently configured.\n']);
    else
        fprintf('COBIRD and MICE are ready: %s\n', mice_root);
    end
end
end
