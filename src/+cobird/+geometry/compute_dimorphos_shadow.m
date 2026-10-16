function shadow_mask = compute_dimorphos_shadow(centroid_Did, dir_Sun_Did, r_Dim_Did, R_Dim)
%COMPUTE_DIMORPHOS_SHADOW Test facet centroids against a spherical-body shadow.
% centroid_Did contains receiving facet centroids [m]; dir_Sun_Did is the
% unit direction toward the Sun. r_Dim_Did and R_Dim give the occluder's
% center [m] and radius [m]. All vectors use the receiving body's frame.
% Either body may receive the shadow. Sunlight is parallel, giving a
% cylindrical shadow; facet incidence is evaluated separately.

to_dimorphos = r_Dim_Did' - centroid_Did;
along_sun = to_dimorphos * dir_Sun_Did;
perpendicular = to_dimorphos - along_sun .* dir_Sun_Did';
shadow_mask = (along_sun > 0) & (vecnorm(perpendicular, 2, 2) < R_Dim);
end
