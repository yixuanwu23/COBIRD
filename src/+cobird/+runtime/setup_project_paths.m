function mice_root = setup_project_paths(project_root, options)
%SETUP_PROJECT_PATHS Add COBIRD and an optional external MICE installation.

arguments
    project_root (1,:) char
    options.require_mice (1,1) logical = true
    options.mice_root (1,:) char = ''
end

project_root = char(java.io.File(project_root).getCanonicalPath());
source_root = fullfile(project_root, 'src');
assert(isfolder(source_root), 'cobird:runtime:setup_project_paths:missingSource', ...
    'COBIRD source directory not found: %s', source_root);

% MATLAB packages are discovered through their parent directory.
addpath(project_root);
addpath(source_root);

configured_root = strtrim(options.mice_root);
if isempty(configured_root)
    configured_root = strtrim(getenv('COBIRD_MICE_ROOT'));
end

if ~isempty(configured_root)
    candidates = {configured_root};
else
    existing_binary = which('mice');
    existing_wrapper = which('cspice_tkvrsn');
    if ~isempty(existing_binary) && ~isempty(existing_wrapper)
        candidates = {fileparts(fileparts(existing_binary))};
    else
        candidates = {};
    end
    if ispc
        platform_directory = 'win64';
    elseif ismac
        platform_directory = 'mac64';
    elseif isunix
        platform_directory = 'linux64';
    else
        error('cobird:runtime:setup_project_paths:unsupportedPlatform', ...
            'This operating system is not supported by the MICE path setup.');
    end
    local_root = fullfile(project_root, 'external', 'mice', platform_directory);
    candidates = [candidates, {local_root, fullfile(local_root, 'mice')}];
end

mice_root = '';
for index = 1:numel(candidates)
    candidate = char(java.io.File(candidates{index}).getCanonicalPath());
    binary = fullfile(candidate, 'lib', ['mice.', mexext]);
    wrapper = fullfile(candidate, 'src', 'mice');
    if isfile(binary) && isfolder(wrapper)
        mice_root = candidate;
        break
    end
end

if isempty(mice_root)
    message = sprintf(['No compatible MICE installation was found. Set ' ...
        'COBIRD_MICE_ROOT to the official MICE root directory. Installation ' ...
        'instructions: https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html']);
    if options.require_mice
        error('cobird:runtime:setup_project_paths:missingMice', '%s', message);
    end
    return
end

% Remove earlier paths below this root before adding one coherent MICE set.
current_paths = strsplit(path, pathsep);
prefix = [mice_root, filesep];
remove_mask = strcmp(current_paths, mice_root) | ...
    strncmp(current_paths, prefix, length(prefix));
if any(remove_mask)
    rmpath(current_paths{remove_mask});
end

addpath(fullfile(mice_root, 'src', 'mice'));
addpath(fullfile(mice_root, 'lib'));
end
