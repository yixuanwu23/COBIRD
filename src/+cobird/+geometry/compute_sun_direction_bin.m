function bin = compute_sun_direction_bin(sun_direction, shadow_cache)
%COMPUTE_SUN_DIRECTION_BIN Map a unit direction to a shadow-cache bin.

azimuth = atan2(sun_direction(2), sun_direction(1));
elevation = asin(max(-1, min(1, sun_direction(3))));
azimuth_bin = max(1, min(shadow_cache.N_az, ...
    floor((azimuth + pi) / shadow_cache.daz) + 1));
elevation_bin = max(1, min(shadow_cache.N_el, ...
    floor((elevation + pi / 2) / shadow_cache.del) + 1));
bin = (elevation_bin - 1) * shadow_cache.N_az + azimuth_bin;
end
