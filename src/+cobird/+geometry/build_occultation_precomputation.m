function occ = build_occultation_precomputation(vertices, faces, normal, no_occ, reabsorb_info)
%BUILD_OCCULTATION_PRECOMPUTATION Precompute data for vectorized radiative recoil.
%
%   occ = cobird.geometry.build_occultation_precomputation(vertices, faces, normal, no_occ, reabsorb_info)
%
%   Returns struct with fields (all in body-fixed frame, constant across time steps):
%     no_occ   : [N×1 logical]  true = no terrain occlusion
%     occ_idx  : [M×1]          global indices of occluded faces
%     x_local  : [M×3]          local x-axis per occluded face
%     y_local  : [M×3]          local y-axis per occluded face
%     Iv_mat   : [M×3]          horizon-integral vector per occluded face

assert(isfield(reabsorb_info, 'method') && strcmp(reabsorb_info.method, 'triangle_horizon'), ...
    'cobird:geometry:build_occultation_precomputation:badReabsorptionMethod', ...
    'Rebuild the reabsorption cache with cobird.geometry.build_reabsorption_cache.');
occ = struct();
occ.no_occ  = no_occ;
occ.occ_idx = find(~no_occ);
M = numel(occ.occ_idx);

if M == 0
    occ.x_local = zeros(0, 3);
    occ.y_local = zeros(0, 3);
    occ.Iv_mat  = zeros(0, 3);
    return
end

idx = occ.occ_idx;

% Local coordinate system (body-fixed, time-invariant)
v1 = vertices(faces(idx, 2), :) - vertices(faces(idx, 1), :);
occ.x_local = v1 ./ vecnorm(v1, 2, 2);
occ.y_local = cross(normal(idx, :), occ.x_local, 2);

% Extract Iv vectors from cell array into [M×3] matrix
occ.Iv_mat = zeros(M, 3);
for k = 1:M
    iv = reabsorb_info.face_occlusion_info{idx(k)};
    occ.Iv_mat(k, :) = iv(:)';
end

end
