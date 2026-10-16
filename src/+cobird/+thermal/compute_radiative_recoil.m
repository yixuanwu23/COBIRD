function [Yarkovsky_force, YORP_moment] = compute_radiative_recoil( ...
    normal, face_centroid, ds, T, params, Q_ref, occ_precomp)
%COMPUTE_RADIATIVE_RECOIL Thermal-emission and reflected-light recoil.
%
%   [F, M] = cobird.thermal.compute_radiative_recoil(normal, centroid, ds, T, params, Q_ref, occ_precomp)
%
%   Force per facet (Sorli et al. 2025, Eq. 6):
%       F_i = -(2/3c) * A_i * n̂_i * (Q_ref_i + ε·σ·T_i⁴)
%   Q_ref is the explicit N_faces-by-1 Lambert-reflected flux [W/m²].
%
%   occ_precomp (required, from cobird.geometry.build_occultation_precomputation):
%     Struct with pre-computed no_occ, occ_idx, x_local, y_local, Iv_mat.
%     When provided, eliminates per-step cellfun + scalar loop overhead.

face_count = size(normal, 1);

if isempty(Q_ref)
    error('cobird:thermal:compute_radiative_recoil:emptyReflectedFlux', ...
        'Q_ref must be an explicit N_faces-by-1 vector.');
end
if ~isequal(size(Q_ref), [face_count, 1])
    error('cobird:thermal:compute_radiative_recoil:reflectedFluxSize', ...
        'Q_ref must be %d-by-1; received %d-by-%d.', ...
        face_count, size(Q_ref, 1), size(Q_ref, 2));
end
if any(~isfinite(Q_ref)) || any(Q_ref < 0)
    error('cobird:thermal:compute_radiative_recoil:invalidReflectedFlux', ...
        'Q_ref must contain finite, nonnegative values.');
end

F      = zeros(face_count, 3);
Q_emit = params.epsilon * params.sigma * T.^4;   % [N×1]
coeff  = -2 / (3 * params.c);

if nargin >= 7 && ~isempty(occ_precomp)
    % Precomputed horizon and reabsorption geometry.
    no_occ  = occ_precomp.no_occ;
    occ_idx = occ_precomp.occ_idx;

    if any(no_occ)
        scale = coeff * ds(no_occ) .* (Q_ref(no_occ) + Q_emit(no_occ));
        F(no_occ, :) = scale .* normal(no_occ, :);
    end

    if ~isempty(occ_idx)
        % Reflected and thermal Lambert photons share this correction.
        Q_out = Q_ref(occ_idx) + Q_emit(occ_idx);
        c_out = -(Q_out .* ds(occ_idx)) / (3 * pi * params.c);
        Iv = occ_precomp.Iv_mat;
        F_out = c_out .* (Iv(:,1) .* occ_precomp.x_local ...
                      + Iv(:,2) .* occ_precomp.y_local ...
                      + Iv(:,3) .* normal(occ_idx, :));

        F(occ_idx, :) = F_out;
    end
else
    error('cobird:thermal:compute_radiative_recoil:noPrecomp', ...
          'occ_precomp is required. Use cobird.geometry.build_occultation_precomputation() to generate it.');
end

M = cross(face_centroid, F, 2);

Yarkovsky_force = sum(F, 1)';
YORP_moment     = sum(M, 1)';

end
