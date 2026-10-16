function hash = compute_file_sha256(path)
%COMPUTE_FILE_SHA256 Return the lowercase SHA-256 digest of one file.

if ~(ischar(path) || (isstring(path) && isscalar(path)))
    error('cobird:io:compute_file_sha256:badPath', ...
        'path must be a character vector or string scalar.');
end
path = char(path);
if ~isfile(path)
    error('cobird:io:compute_file_sha256:fileNotFound', 'File not found: %s', path);
end

fid = fopen(path, 'rb');
if fid < 0
    error('cobird:io:compute_file_sha256:openFailed', 'Unable to read file: %s', path);
end
cleaner = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');

digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(bytes, 'int8'));
hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end
