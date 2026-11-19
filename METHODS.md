# Per_9(0): how the model was computed

## The curve

Per_9(0) is the curve of quadratic rational maps with a marked critical point of period 9. In the chart

    f(z) = (z - x)(z - w)/z^2,   w = (x + y - 1)/(x - 1),

the marked cycle is 0 -> inf -> 1 -> y -> a_4 -> ... -> a_7 -> x -> 0 (y = a_3, x = a_8), and the curve is P(x, y) = 0
(`per9_xy_poly.txt`), of degree 42 in x. The free critical point is omega = 2xw/(x + w). The compact curve has genus 188,
252 polynomial centers, 156 punctures (144 with m = 1 and 12 with m = 2 in the census) and at least 6 real circles.

## Steps

1. **Branched cover** (`branch.json`, computed in an earlier probe with the same root tracker `trackmp.py`): y : Per_9(0) -> sphere
   has degree 42; 419 candidate branch points, 408 with non-trivial monodromy, total ramification 458; Riemann-Hurwitz gives genus 188.
2. **Centers** (`make_centres.py`) and the explorer data (`explorer/build_data9.py` -> `explorer/data.json`): the 252 centers located on the sheets.
3. **Mesh of the y-sphere** (`surface3d/make_sphere_mesh.py`): 96 817 vertices, 193 630 triangles, 408 branch vertices; radii near the nodes in `nodes.json`.
4. **Sheets and gluing** (`surface3d/glue.py`, which uses `cutglue.py`, then `surface3d/assemble.py`): cut system on the sphere, the 42 sheets glued into one closed
   surface of 4 065 856 vertices and 8 132 460 triangles, Euler characteristic -374.
5. **Uniformization** (`surface3d/uniformize.py`, then `surface3d/refine.py` to 160-bit accuracy): discrete conformal change of metric (Bobenko-Pinkall-Springborn)
   to a hyperbolic metric with angle sum 2 pi at every vertex.
6. **Development and fundamental polygon** (`surface3d/develop.py`, with `fuchs.py`): development into the Poincare disk and the Dirichlet polygon about a base point,
   2250 = 12g - 6 sides.
7. **Values on the mesh** (`surface3d/fold_vertices9.py`, then `vertex_x.py`, `vertex_x_polish.py`, `vertex_x_finish.py`, `vertex_x_inner.py`): every mesh vertex folded into
   the polygon, and the coordinate x on every sheet at every vertex.
8. **Charts for the viewer** (`surface3d/chart_export9.py`, then `surface3d/pack_grids9.py`): 804 charts storing the two zeros (w, x) of the map (y = (1 - x)(1 - w)),
   96 x 96 grids; points with |y| < 0.04 are masked. The grids are packed losslessly into `chart_grids_0..4.dxs` (format in `pack_grids9.py`).
9. **Polygon and markers for the viewer** (`surface3d/conformal_export9.py` -> `conformal_data.js`; `surface3d/marker_info9.py` -> `marker_info.json` for the centers).
10. **Punctures** (`puncture_place.py` -> `puncture_place.json`; `puncture_markers.py` adds them to `explorer/data.json`, `conformal_data.js`, `marker_info.json`): each
    puncture is placed as the boundary point of a hyperbolic component reached along an internal ray.
11. **Overlays** (`surface3d/g9_edges.py`: the graph G_9; `surface3d/trace_ovals.py`, `surface3d/oval_check9.py`, `surface3d/oval_export9.py`: the real circles that pass
    the check) -> `overlay_data.json`.

## Checks

Each check was written down before it was run. A coarser mesh (not included) was computed alongside for the convergence checks.

| check | result |
|---|---|
| Riemann-Hurwitz from the monodromy | genus 188 (the predicted value) |
| mesh closed and orientable, Euler characteristic -374 | chi = -374 on both meshes; 0 closure failures |
| the 252 centers lie on the curve | max relative residual 3.6e-17 |
| angle sum 2 pi at every vertex | 3.8e-12 in double precision; 1.5e-34 after refinement |
| total hyperbolic area = 4 pi (g - 1) = 748 pi | difference 1.9e-9 |
| side pairings are isometries; vertex cycles multiply to the identity | side residual 3.4e-25, cycles 1.5e-16 |
| Dirichlet polygon has area 748 pi and its sides pair up | 2250 sides; difference -1.2e-13; every partner is the inverse element |
| mesh convergence: the 31 626 distances between the 252 centers, coarse vs fine mesh | with each center inserted as a mesh vertex: relative difference max 3.0e-2, median 5.6e-4 (passes; threshold 5%) |
| the punctures are distinct places | see the limitations: the first method placed 121 of 156; the placement used (internal rays) places all 156, all distinct |
| the real circles are closed geodesics (axis distance < 2e-2, length ratio within 1e-2, both shrinking on the fine mesh) | circles 0, 3, 4, 5 pass (axis distance 7.3e-3, 1.3e-2, 4.6e-3, 4.4e-3); circles 1 and 2 fail (5.9e-2, 1.8e-2) and are not drawn |
| charts reproduce the mesh values | chordal error of y over 3991 evaluations: median 5.0e-6, 99th percentile 1.2e-4, max 8.3e-4; fraction above 1e-3: 0 |

## Known limitations

- **Mesh convergence**: with the centers placed by barycentric coordinates in their mesh triangles, one pair (centers 135 and 136, on different real branches next to a
  cluster of branch points near y = 1) differs by 8.0% between the meshes, above the 5% threshold. With the centers inserted as mesh vertices the check passes; that
  placement is the one recorded.
- **Punctures**: the first placement method (continuing along the sheets) located 156 punctures as only 127 distinct places, with 35 continuation failures. The
  punctures shown are placed by internal rays instead: 156 of 156, all distinct; where both methods placed a puncture (119 cases) they agree.
- **Real circles**: circles 1 and 2 fail the geodesic test; their largest deviations are near y about 0.60 (circle 1) and y about 0.99 (circle 2). They are not
  drawn. The number of real circles is a lower bound: 8 real arcs of the census did not land.
- **Charts near y = 0**: points with |y| < 0.04 are masked (shown in a neutral colour); there the sheets nearly coincide and single-precision charts cannot separate them.
