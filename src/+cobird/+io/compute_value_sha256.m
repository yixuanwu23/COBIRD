function hash = compute_value_sha256(value)
%COMPUTE_VALUE_SHA256 Fingerprint numeric configuration and geometry data.
% Field order is ignored. Array class, size, and sparse indices are retained.
digest = java.security.MessageDigest.getInstance('SHA-256');
local_update(digest, value);
hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end

function local_update(digest, value)
header = unicode2native(jsonencode({class(value), size(value)}), 'UTF-8');
digest.update(typecast(uint8(header), 'int8'));
if isstruct(value)
    names = sort(fieldnames(value));
    local_update(digest, names);
    for index = 1:numel(value)
        for field = 1:numel(names)
            local_update(digest, value(index).(names{field}));
        end
    end
elseif iscell(value)
    for index = 1:numel(value)
        local_update(digest, value{index});
    end
elseif ischar(value) || isstring(value)
    bytes = unicode2native(jsonencode(value), 'UTF-8');
    digest.update(typecast(uint8(bytes), 'int8'));
elseif isnumeric(value) || islogical(value)
    assert(isreal(value), 'cobird:io:compute_value_sha256:complex', ...
        'Complex configuration arrays are not supported.');
    digest.update(int8(issparse(value)));
    if issparse(value)
        [row, column, entries] = find(value);
        local_update(digest, row);
        local_update(digest, column);
        local_update(digest, entries);
    elseif ~isempty(value)
        if islogical(value)
            value = uint8(value);
        end
        digest.update(typecast(value(:), 'int8'));
    end
else
    error('cobird:io:compute_value_sha256:type', ...
        'Unsupported configuration type: %s.', class(value));
end
end
