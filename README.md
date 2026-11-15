# Per_6(0): the computation behind the model

Scripts and data that produce the uniformized model of Per_6(0) shown in the
[Per_n(0) atlas](https://callieshlb.github.io/per_n-viewer/). This branch is not part of the website (the website is the
`main` branch). `METHODS.md` describes the steps, the checks and their results, and the known limitations.

## Layout

- `per6/` — the computation. Top level: branched cover (branch points, sheets, monodromy). `explorer/`: export of the curve's
  coefficients and census data. `surface3d/`: mesh, uniformization, development into the disk, charts, overlays.
  `surface3d/coarse/`: scripts for the coarser mesh used in the convergence check (its data is not included).
- `inputs/` — small files the scripts read that were computed earlier: the equation of the curve (`per6_xy_poly.txt`), the census
  of polynomial centers, punctures and real-circle arcs with the helper module that reads it (`E90-real-circles/`), and the
  graph G_6 (`E91-Gn-graph/`).

## Order of the scripts

1. `per6/atlas_branch.py` (uses `atlas_lib.py`) -> `branch.json`; `per6/atlas_render.py` -> `atlas_centers.json`
2. `per6/explorer/export_data.py` -> `explorer/data.json`
3. `per6/surface3d/make_sphere_mesh.py` -> `sphere_mesh.json`; `edge_perms.js` (node) -> `edge_perms.json`; `assemble.py` -> `surface_mesh.json`
4. `uniformize.py` -> `uniform_u.npy`; `develop.py` -> `develop_pre.npz`, `develop.npz`; `dirichlet_only.py` reruns only the Dirichlet polygon
5. `colour_export.py` -> `surface_data.js`; `conformal_export.py` -> `conformal_data.js`; `chart_export.py` -> `chart_grids.bin`, `chart_meta.json`, `chart_map.png`
6. `oval_check.py` -> `oval_points.json`, `oval_lifts.npz`, `oval_check.json`; `overlay_export.py` -> `overlay_data.json`

`meshloc.py` is a shared helper (locates a point of the curve in the mesh and in the disk).

## Data

| file | size | what |
|---|---|---|
| `per6/branch.json`, `per6/atlas_centers.json`, `per6/explorer/data.json` | small | branch points and monodromy; centers and punctures on the sheets; coefficients |
| `surface3d/sphere_mesh.json` | 10.5 MB | mesh of the y-sphere |
| `surface3d/edge_perms.json` | 24.1 MB | roots at each vertex and sheet permutations along mesh edges |
| `surface3d/surface_mesh.json` | 23.2 MB | the glued surface mesh |
| `surface3d/uniform_u.npy` | 3.0 MB | conformal factors |
| `surface3d/develop_pre.npz`, `develop.npz` | 36.3 MB | development into the disk, group generators, Dirichlet polygon |
| `surface3d/surface_data.js`, `tiling_data.js`, `conformal_data.js` | 58.3 MB | per-vertex exports; the last is what the viewer loads |
| `surface3d/chart_grids.bin`, `chart_meta.json`, `chart_map.png` | 10.6 MB | charts (what the viewer loads) |
| `surface3d/overlay_data.json`, `oval_points.json`, `oval_lifts.npz`, `oval_check.json` | small | graph edges and real circles |

## Notes

- `conformal_export.py` reads the per-vertex coordinates from `tiling_data.js`. That file is included, but the script that wrote
  it is not: it belongs to a separate construction (a pants decomposition of the surface) that is not part of this release.
- The scripts originally used absolute paths on the machine they ran on. These were replaced by paths relative to each script.
  After that change only the first stage was tried (loading the curve and branch data); the full pipeline was not re-run.
- Requirements: Python with numpy, scipy, sympy, mpmath, Pillow; node for `edge_perms.js`.
- In the file names and script headers, E90, E91, E92 are internal run numbers: E90 = the census, E91 = the graph G_6, E92 = this computation.

Working files of ongoing research. No licence has been chosen yet; please ask before reusing.
