# Per_7(0): the computation behind the model

Scripts and data that produce the uniformized model of Per_7(0) shown in the
[Per_n(0) atlas](https://callieshlb.github.io/per_n-viewer/). This branch is not part of the website (the website is the
`main` branch). `METHODS.md` describes the steps, the checks and their results, and the known limitations.

## Layout

- `per7/` — branched cover: `disc.py` (candidate branch points), `trackmp.py` (256-bit root continuation), `atlas_branch.py`
  (monodromy, genus) -> `branch.json`; `atlas_points.py` -> `atlas_points.json`, `explorer/data.json` (centers and punctures on the sheets).
- `per7/explorer/` — `per7_core.js` (double-precision dynamics and root tracker used by the viewer), `test_core.js` (checks that every center is a polynomial center), `data.json`.
- `per7/surface3d/` — mesh, sheets, uniformization, development, charts, overlays: scripts and their outputs.
- `inputs/` — small files computed earlier: the equation of the curve and the census of centers, punctures and real-circle arcs
  with its helper module (`E90-real-circles/`), and the graph G_7 (`E91-Gn-graph/`).

## Order (run in `per7/surface3d/` unless noted)

1. `per7/disc.py`; `per7/atlas_branch.py` -> `branch.json`; `per7/atlas_points.py` -> `atlas_points.json`, `explorer/data.json`
2. `make_sphere_mesh.py` -> `sphere_mesh.json`
3. `edge_perms.py` -> `edge_perms.json`; `assemble.py .` -> `surface_mesh.json`
4. `uniformize.py` -> `uniform_u.npy`; `refine.py` -> `uniform_u_hi.npy`, `uniform_u_lo.npy`
5. `develop.py mp2` -> `develop_pre.npz`; `dirichlet_only.py` -> `develop.npz`
6. `trace_ovals.py` -> `oval_points.json`; `oval_check.py` -> `oval_lifts.npz`, `oval_check.json`
7. `colour_export.py` -> `surface_data.json`; `conformal_export.py` -> `conformal_data.js`
8. `chart_export.py` -> `chart_grids.bin`, `chart_meta.json`, `chart_map.png`; `chart_map_fill.py` fills the edge pixels of `chart_map.png`
9. `overlay_export.py` -> `overlay_data.json`; `marker_info.py` -> `marker_info.json`

`meshloc.py` (locating a point of the curve in the mesh and the disk) and `fuchs.py` (multiprecision group and Dirichlet
polygon code) are shared helpers. The files the site loads are `conformal_data.js`, `chart_grids.bin`, `chart_meta.json`,
`chart_map.png`, `overlay_data.json`, `marker_info.json` and `explorer/data.json`.

## Notes

- Not included: the coarse mesh used for the convergence checks (`make_sphere_mesh.py coarse` builds its sphere mesh; the later scripts take the mesh directory),
  and the convergence script itself.
- E90, E91, E93 in names are internal run numbers: E90 = the census, E91 = the graph G_7, E93 = this computation.
- Nine paths were changed so the scripts find `inputs/` in this layout. The pipeline was not re-run after that change.
- Requirements: Python with numpy, scipy, sympy, mpmath, gmpy2, Pillow; node for `explorer/test_core.js`.

Working files of ongoing research. No licence has been chosen yet; please ask before reusing.
