# Per_6(0): the computation behind the model

Scripts and data that produce the uniformized model of Per_6(0) shown in the
[Per_n(0) atlas](https://callieshlb.github.io/per_n-viewer/). This branch is not part of the website (the website is the
`main` branch). `METHODS.md` describes the steps, the checks and their results, and the known limitations.

**Version 2** (2026-11-18): mesh refined at the punctures, extra chart at the cone point 6/58. This is the model the site shows.
Version 1 (the first mesh) is the previous commit of this branch.

## Layout

- `per6/` — branched cover (branch points, sheets, monodromy): `atlas_lib.py`, `atlas_branch.py`, `atlas_render.py`, with `branch.json`, `atlas_centers.json`.
- `per6/explorer/` — the curve's coefficients (`export_data.py` -> `data.json`) and the root tracker used to label sheets (`per6_core.js`; `per6_core_mesh.js` is the same with a smaller minimum step, for the mesh).
- `per6/surface3d/` — mesh, sheets, uniformization, development, charts, overlays: scripts and their outputs.
- `per6/cone/` — the extra chart at 6/58. `chart_grids.bin` and `chart_meta.json` here (25 charts) are the files the site loads; those in `surface3d/` (24 charts) are the step before.
- `inputs/` — small files computed earlier that the scripts read: the equation of the curve, the census of centers, punctures and real-circle arcs with its helper module (`E90-real-circles/`), and the graph G_6 (`E91-Gn-graph/`).

## Order (run in `per6/surface3d/` unless noted)

1. `per6/atlas_branch.py` -> `branch.json`; `per6/atlas_render.py` -> `atlas_centers.json`; `per6/explorer/export_data.py` -> `explorer/data.json`
2. `make_sphere_mesh.py` -> `sphere_mesh.json`
3. `node edge_perms.js` -> `edge_perms.json`; then `edge_perms_mp.py 4` repairs it in 40 digits near y = 0, 1/2, 1
4. `assemble.py` -> `surface_mesh.json`
5. `uniformize.py` -> `uniform_u.npy`
6. `develop.py mp` -> `develop_pre.npz` (uses `devmp.py`); `dirichlet_only.py` -> `develop.npz`
7. `oval_check.py` -> `oval_points.json`, `oval_lifts.npz`, `oval_check.json`
8. `colour_export.py` -> `surface_data.js`; `vertex_fields.py` -> `vertex_fields.js`
9. `chart_export.py` -> `chart_grids.bin`, `chart_meta.json`, `chart_map.png`; `chart_map_fill.py` fills the edge pixels of `chart_map.png`
10. `conformal_export.py` -> `conformal_data.js`; `overlay_export.py` -> `overlay_data.json`
11. cone at 6/58: `python3 ../cone/chart_export64.py` (from `surface3d/`) -> `cone/chart_grids64.npy`; then `cone/cone_data.py` -> the 25-chart files

`meshloc.py` is a shared helper (locates a point of the curve in the mesh and in the disk).

## Notes

- Not included: `cone/chart_grids64.npy` (27 MB, an intermediate of step 11), and the inputs of `vertex_fields.py`, which
  copies a click-info label ("piece") from the version-1 mesh; the output `vertex_fields.js` is included.
- Script comments refer to `REGISTRATION.md`, the internal record of the computation; it is not included (`METHODS.md` summarises it).
  E90, E91, E92, E98, E99 in names are internal run numbers: E90 = the census, E91 = the graph G_6, E92 = the version-1 model, E99 = the refinement.
- Three paths were changed so the scripts find `inputs/` and `surface3d/` in this layout. The pipeline was not re-run after that change.
- Requirements: Python with numpy, scipy, sympy, mpmath, Pillow; node for `edge_perms.js` and `edge_perms_mp.py` (which calls node).

Working files of ongoing research. No licence has been chosen yet; please ask before reusing.
