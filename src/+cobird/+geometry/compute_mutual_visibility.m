function [vis_i, vis_j, cos_i, cos_j, dist_ij] = compute_mutual_visibility( ...
    centroid_Did, normal_Did, reabsorb_Did, ...
    centroid_Dim_inDid, normal_Dim_inDid, reabsorb_Dim, ...
    r_Dim_inDid, R_did, R_dim, sh_Did, sh_Dim, dir_to_Did_inDimBody)
%COMPUTE_MUTUAL_VISIBILITY Unified centerline mutual-radiation geometry.
% Both normal arrays are expressed in the primary-body frame.
% d_hat points from the Didymos center to the Dimorphos center. The emitting
% and receiving cosines are n_Did*d_hat and n_Dim*(-d_hat), respectively.
% Both shadow tables are queried along this same line in their own frames.
% The caller supplies dir_to_Did_inDimBody = R_Dim_to_Did'*(-d_hat).
%
% All strictly positive cosines are retained. Surviving facets form
% separate lists for the primary and secondary:
%   q_Dim(j) = epsilon*cos_j(j) * ...
%       sum(cos_i .* area_Did(vis_i) .* exitance_Did(vis_i)) / (pi*d^2).
% The reverse exchange uses the same geometry. Cost is O(N_Did+N_Dim).
% Pair distances and directions use the common centerline approximation;
% visibility is evaluated along that direction on each body.
%
% Centroid and radius arguments are unused in this calculation.
% All outputs are empty when either body has no surviving facets.

r_vec = r_Dim_inDid(:);
dist_COM = norm(r_vec);
assert(numel(r_vec) == 3 && all(isfinite(r_vec)) && ...
    isfinite(dist_COM) && dist_COM > 0, ...
    'cobird:geometry:compute_mutual_visibility:invalidSeparation', ...
    'A finite, nonzero three-component center separation is required.');
dir_Did = r_vec / dist_COM;
dir_Dim = dir_to_Did_inDimBody(:)';
dir_Dim_norm = norm(dir_Dim);
assert(numel(dir_Dim) == 3 && all(isfinite(dir_Dim)) && ...
    isfinite(dir_Dim_norm) && dir_Dim_norm > 0, ...
    'cobird:geometry:compute_mutual_visibility:invalidDimorphosDirection', ...
    'A finite, nonzero Didymos direction in the Dimorphos frame is required.');
dir_Dim = dir_Dim / dir_Dim_norm;

cos_Did = normal_Did * dir_Did;
cos_Dim = normal_Dim_inDid * (-dir_Did);
cand_Did = find(cos_Did > 0);
cand_Dim = find(cos_Dim > 0);
if isempty(cand_Did) || isempty(cand_Dim)
    vis_i=[]; vis_j=[]; cos_i=[]; cos_j=[]; dist_ij=[];
    return
end

surv_Did = local_shadow_horizon_filter(cand_Did, reabsorb_Did, sh_Did, dir_Did);
surv_Dim = local_shadow_horizon_filter(cand_Dim, reabsorb_Dim, sh_Dim, dir_Dim);
if isempty(surv_Did) || isempty(surv_Dim)
    vis_i=[]; vis_j=[]; cos_i=[]; cos_j=[]; dist_ij=[];
    return
end
vis_i = uint32(surv_Did(:));
vis_j = uint32(surv_Dim(:));
cos_i = cos_Did(surv_Did);
cos_j = cos_Dim(surv_Dim);
dist_ij = dist_COM;
end

function surv = local_shadow_horizon_filter(cand, reabsorb, sh, dir_body)
% Direction lookup is constant-time per facet, O(number of candidates) total.
az = atan2(dir_body(2), dir_body(1));
el = asin(max(-1, min(1, dir_body(3))));
az_bin = max(1, min(sh.N_az, floor((az + pi) / sh.daz) + 1));
el_bin = max(1, min(sh.N_el, floor((el + pi/2) / sh.del) + 1));
bin = (el_bin - 1) * sh.N_az + az_bin;
has_cand = reabsorb.has_candidates(cand);
cand_with_tri = cand(has_cand);
cand_no_tri = cand(~has_cand);
if isempty(cand_with_tri)
    surv = cand_no_tri;
    return
end
row_idx = sh.occ_rev_map(cand_with_tri);
assert(all(row_idx > 0), 'cobird:geometry:compute_mutual_visibility:invalidShadowMap', ...
    'The shadow-cache occ_rev_map is inconsistent with reabsorption metadata.');
blocked = sh.shadow_table(row_idx, bin);
surv = [cand_no_tri; cand_with_tri(~blocked)];
end
