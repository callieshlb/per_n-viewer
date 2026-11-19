# Per_8(0): the computation behind the model

Scripts and data that produce the uniformized model of Per_8(0) shown in the
[Per_n(0) atlas](https://callieshlb.github.io/per_n-viewer/). This branch is not part of the website (the website is the
`main` branch). `METHODS.md` describes the steps, the checks and their results, and the known limitations.

## Layout

- `per8/` — the curve (`per8_xy_poly.txt`) and the branched cover: `disc.py`, `trackmp.py` (256-bit root continuation),
  `atlas_branch.py` -> `branch.json`; `nodes.json` (mesh radii at the nodes); `atlas_points.py` -> `atlas_points.json`, `explorer/data.json`.
- `per8/explorer/` — `per8_core.js` (double-precision dynamics and root tracker used by the viewer), `test_core.js`, `data.json`.
- `per8/surface3d/` — mesh, sheets, uniformization, development, charts, overlays: scripts and their outputs.
- `inputs/` — files computed earlier: the census of centers, punctures and real-circle arcs with its helper module
  (`E90-real-circles/`), and the graph G_8 (`E91-Gn-graph/`).

## Order (run in `per8/surface3d/` unless noted)

1. `per8/disc.py`; `per8/atlas_branch.py` -> `branch.json`; `per8/atlas_points.py` -> `atlas_points.json`, `explorer/data.json`
2. `make_sphere_mesh.py` -> `sphere_mesh.json`
3. `edge_perms.py` -> `edge_perms.json`; `assemble.py .` -> `surface_mesh.json`
4. `uniformize.py` -> `uniform_u.npy`; `refine.py` -> `uniform_u_hi.npy`, `uniform_u_lo.npy`
5. `develop.py mp2` -> `develop_pre.npz`; `u4_bfs.py` then `dirichlet_only.py` -> `develop.npz`
6. `trace_ovals.py` -> `oval_points.json`; `oval_check.py` -> `oval_lifts.npz`, `oval_check.json`
7. `colour_export.py` -> `surface_data.json`; `conformal_export.py` -> `conformal_data.js`
8. `chart_export.py` -> `chart_grids_0..2.bin`, `chart_meta.json`, `chart_map.png`; `chart_place_all.py` (all placements, into `chart_meta.json`);
   `chart_polish.py` (grid nodes onto the curve, rewrites the grids); `chart_map_fill.py` (edge pixels of `chart_map.png`)
9. `overlay_export.py` -> `overlay_data.json`; `marker_info.py` -> `marker_info.json`

`meshloc.py` (locating a point of the curve in the mesh and the disk) and `fuchs.py` (multiprecision group and Dirichlet
polygon code) are shared helpers. The files the site loads are `conformal_data.js`, `chart_grids_0..2.bin`, `chart_meta.json`,
`chart_map.png`, `overlay_data.json`, `marker_info.json` and `explorer/data.json`.

## Notes

- Not included: the coarse mesh used for the convergence checks (`make_sphere_mesh.py coarse` builds its sphere mesh; the
  later scripts take the mesh directory) and the convergence script itself.
- E90, E91, E94 in names are internal run numbers: E90 = the census, E91 = the graph G_8, E94 = this computation.
- Four paths were changed so the scripts find `inputs/` in this layout. The pipeline was not re-run after that change.
- Requirements: Python with numpy, scipy, sympy, mpmath, gmpy2, Pillow; node for `explorer/test_core.js`.
- Three files are larger than 50 MB (the two largest chart grid parts and `surface_mesh.json`); all are below GitHub's 100 MB limit.

Working files of ongoing research. No licence has been chosen yet; please ask before reusing.
