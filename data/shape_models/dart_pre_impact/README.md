# Official DART pre-impact shape models

This directory contains unmodified NASA Planetary Data System products used
by the COBIRD pre-impact Didymos example:

| Body | PDS product | Status | Vertices | Facets |
|---|---|---|---:|---:|
| Didymos | `didymos_g_9309mm_spc_obj_0000n00000_v003` | v003 final | 24,578 | 49,152 |
| Dimorphos | `dimorphos_g_1940mm_spc_obj_0000n00000_v004` | v004 final pre-impact | 24,578 | 49,152 |

Each OBJ is paired with its original PDS4 XML label. The files are retained
byte-for-byte; COBIRD converts vertex coordinates from kilometres to metres
only after loading them into memory. The source products are not recentered,
rotated, rescaled, decimated, or rewritten.

## Geometry caches

`cobird.geometry.prepare_model_cache` expects an OBJ already expressed in
metres. Do not pass these kilometre-based PDS files directly to that entry
point. The lower-level builders accept in-memory metre-based vertices.
The following prepares both bodies with terrain caches in a separate output
directory. Check the source files against `SHA256SUMS` before preprocessing.
Building full-mesh terrain caches takes substantially more work than the
bundled quick example.

```matlab
project_root = startup_cobird;
data_root = fullfile(project_root, "data", "shape_models", "dart_pre_impact");
cache_root = fullfile(project_root, "local_cache");
if ~isfolder(cache_root), mkdir(cache_root); end
model_names = { ...
    'didymos_g_9309mm_spc_obj_0000n00000_v003', ...
    'dimorphos_g_1940mm_spc_obj_0000n00000_v004'};
bodies = cell(1,2);
for index = 1:2
    name = model_names{index};
    [vertices_km, faces] = cobird.geometry.load_obj_mesh( ...
        fullfile(data_root, [name, '.obj']));
    vertices = vertices_km * 1000;
    [normal, area, centroid] = cobird.geometry.compute_face_geometry(vertices, faces);

    cobird.geometry.build_reabsorption_cache(vertices, faces, cache_root, name);
    loaded = load(fullfile(cache_root, [name, '_reabsorb_info.mat']), 'reabsorb_info');
    reabsorb = loaded.reabsorb_info;
    cobird.geometry.build_shadow_horizon_cache( ...
        vertices, faces, normal, centroid, reabsorb, 360, 180, cache_root, name);
    shadow = cobird.geometry.load_shadow_horizon_cache(cache_root, name, 360, 180);
    no_occ = cellfun(@isempty, reabsorb.face_occlusion_info);
    occultation = cobird.geometry.build_occultation_precomputation( ...
        vertices, faces, normal, no_occ, reabsorb);

    vf_parameters = struct('F_threshold', 1e-4, 'ray_offset_m', 1e-6);
    cobird.geometry.build_view_factor_cache( ...
        vertices, faces, normal, area, centroid, cache_root, name, vf_parameters);
    loaded = load(fullfile(cache_root, [name, '_VF.mat']), 'VF');
    bodies{index} = struct('vertices', vertices, 'faces', faces, ...
        'normal', normal, 'area', area, 'centroid', centroid, ...
        'radius_m', max(vecnorm(vertices, 2, 2)), 'reabsorb', reabsorb, ...
        'shadow', shadow, 'occultation', occultation, 'view_factor', loaded.VF);
end
primary = bodies{1};
secondary = bodies{2};
```

Use these bodies with `cobird.prepare_binary(..., terrain_heating=true)` and
the physical parameters, initial state, and environment for the chosen case.
See `help cobird.prepare_binary` for those inputs. Each builder writes its
named cache products; keep caches from different geometries or settings in
separate directories. The original OBJ/XML files remain unchanged.

## Provenance and citation

- DART Shapemodel Archive Bundle, DOI
  [10.26007/96fn-p578](https://doi.org/10.26007/96fn-p578).
- Didymos v003 collection, DOI
  [10.26007/bm57-x327](https://doi.org/10.26007/bm57-x327).
- Dimorphos v004 collection, DOI
  [10.26007/0nss-vd15](https://doi.org/10.26007/0nss-vd15).

The archive citation lists T. Daly, O. Barnouin, C. Ernst, H. Nair,
R. Espiritu, and D. Waller (2023). Cite the PDS bundle or the most specific
collection used, in addition to citing COBIRD. File integrity values are
recorded in `SHA256SUMS`.

These PDS data products are external scientific data and are not covered by
the COBIRD software license.
