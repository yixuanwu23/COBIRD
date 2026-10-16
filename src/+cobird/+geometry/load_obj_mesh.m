function [vertices, faces] = load_obj_mesh(filename)
%LOAD_OBJ_MESH Read a triangular OBJ mesh.

if ~isfile(filename)
    error('cobird:geometry:load_obj_mesh:fileOpenFailed', ...
        'Unable to open OBJ file: %s', filename);
end

lines = strtrim(regexprep(readlines(filename), '#.*$', ''));
lines = regexprep(lines, '\s+', ' ');
vertex_lines = lines(startsWith(lines, "v "));
face_lines = lines(startsWith(lines, "f "));
assert(~isempty(vertex_lines) && ~isempty(face_lines), ...
    'cobird:geometry:load_obj_mesh:emptyMesh', ...
    'OBJ file contains no vertices or faces: %s', filename);

vertices = zeros(numel(vertex_lines), 3);
for vertex_index = 1:numel(vertex_lines)
    vertex_fields = split(vertex_lines(vertex_index));
    assert(numel(vertex_fields) >= 4, ...
        'cobird:geometry:load_obj_mesh:invalidVertex', ...
        'OBJ file contains an invalid vertex record: %s', filename);
    vertices(vertex_index, :) = str2double(vertex_fields(2:4)).';
end

face_fields = split(face_lines);
if isscalar(face_lines)
    face_fields = face_fields.';
end
assert(size(face_fields, 2) == 4, ...
    'cobird:geometry:load_obj_mesh:nonTriangularFace', ...
    'Only triangular OBJ faces are supported: %s', filename);
face_tokens = face_fields(:, 2:4);
faces = str2double(extractBefore(face_tokens + "/", "/"));

assert(all(isfinite(vertices), 'all') && all(isfinite(faces), 'all'), ...
    'cobird:geometry:load_obj_mesh:invalidNumber', ...
    'OBJ file contains invalid numeric values: %s', filename);
assert(all(faces >= 1, 'all') && all(faces == fix(faces), 'all'), ...
    'cobird:geometry:load_obj_mesh:invalidFaceIndex', ...
    'OBJ file contains invalid face indices: %s', filename);
assert(all(faces <= size(vertices, 1), 'all'), ...
    'cobird:geometry:load_obj_mesh:faceIndexOutOfRange', ...
    'OBJ file contains a face index outside the vertex table: %s', filename);
end
