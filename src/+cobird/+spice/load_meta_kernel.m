function load_meta_kernel(meta_kernel)
%LOAD_META_KERNEL Load relative kernel paths from a meta-kernel directory.

assert(isfile(meta_kernel), 'cobird:spice:load_meta_kernel:fileNotFound', ...
    'SPICE meta-kernel not found: %s', meta_kernel);
[kernel_dir, kernel_name, kernel_ext] = fileparts(meta_kernel);
original_dir = pwd;
restore_dir = onCleanup(@() cd(original_dir));
cd(kernel_dir);
cspice_furnsh([kernel_name, kernel_ext]);
clear restore_dir
end
