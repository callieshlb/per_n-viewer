# Per_8(0): how the model was computed

## The curve

Per_8(0) is the curve of quadratic rational maps with a marked critical point of period 8. In the chart

    f(z) = (z - x)(z - w)/z^2,   w = (x + y - 1)/(x - 1),

the marked cycle is 0 -> inf -> 1 -> y -> a_4 -> a_5 -> a_6 -> x -> 0 (y = a_3, x = a_7), and the curve is the main factor
P(x, y) = 0 of the closing equation, with deg_x P = deg_y P = 20 (`per8_xy_poly.txt`). The free critical point is
omega = 2xw/(x + w). The compact curve has genus 69, 120 polynomial centers, 74 punctures and 5 real circles.

## Steps

1. **Branched cover.** y : Per_8(0) -> sphere has degree 20. Candidate branch points are the roots of the discriminant of P in x
   and of its leading coefficient (`disc.py`). Roots are continued in 256-bit arithmetic (`trackmp.py`), because they cluster
   strongly near the nodes. Monodromy of a loop around each candidate (`atlas_branch.py`): total ramification 176, nodes
   (trivial monodromy, colliding roots) at y = 0 and at the two roots of y^2 + y - 1, no pole points.
   Riemann-Hurwitz: g = (176 - 38)/2 = 69.
2. **Centers and punctures on the sheets** (`atlas_points.py`): the 120 centers and 74 puncture loop points of the census.
3. **Mesh of the y-sphere.** A refined icosphere with branch points as vertices, rings around branch points, finer grids where
   they cluster, and no vertices within 0.08 of y = 0 and 0.01 of the other two nodes (`surface3d/make_sphere_mesh.py`, radii in `nodes.json`).
4. **Sheets and gluing.** Roots at every vertex and the sheet permutation along every edge, by 256-bit continuation
   (`surface3d/edge_perms.py`); the sheets glued into one closed surface (`surface3d/assemble.py`): 844 924 vertices,
   1 690 120 triangles, Euler characteristic -136.
5. **Uniformization.** Discrete conformal change of metric (Bobenko-Pinkall-Springborn) to a hyperbolic metric with angle sum
   2 pi at every vertex (`surface3d/uniformize.py`), refined to 160-bit accuracy (`surface3d/refine.py`).
6. **Development and fundamental polygon.** Development into the Poincare disk in 160 bits (`surface3d/develop.py mp2`,
   `devmp.py`, `devmp2.py`). The Dirichlet polygon (822 = 12g - 6 sides, 274 vertex cycles) is computed with all geometry in
   200 bits, seeded by a nearest-lift walk over the mesh (`surface3d/u4_bfs.py`, `surface3d/dirichlet_only.py`, `surface3d/fuchs.py`).
7. **Charts for the viewer.** (x, y) as functions on the disk on 318 overlapping charts, fitted to the mesh by local holomorphic
   least squares (`surface3d/chart_export.py`). Every placement of a chart covering the polygon is exported (12 499;
   `surface3d/chart_place_all.py`), so the viewer can choose, per view, the chart each pixel is deepest in. Every grid node is
   then moved onto the curve: y kept, x replaced by the nearest root of P(., y) (`surface3d/chart_polish.py`). Edge pixels of the
   chart map are assigned to the nearest chart (`surface3d/chart_map_fill.py`).
8. **Overlays.** Real circles (`surface3d/trace_ovals.py`, `surface3d/oval_check.py`), the graph G_8, centers and punctures
   (`surface3d/overlay_export.py`, `surface3d/conformal_export.py`, `surface3d/marker_info.py`).

## Checks

Each check was written down before it was run. A coarser mesh (not included) was computed alongside for the convergence checks.

| check | result |
|---|---|
| Riemann-Hurwitz from the monodromy; product of the loops = identity; monodromy group transitive | genus 69; both hold |
| mesh closed and orientable, Euler characteristic -136 (second route to the genus) | chi = -136 on both meshes |
| the 120 centers lie on the curve | max relative residual 1.5e-17 |
| the 74 punctures are 74 distinct places | passes; 0 continuation failures |
| angle sum 2 pi at every vertex | 3.7e-13 |
| total hyperbolic area = 4 pi (g - 1) = 272 pi | difference 6.8e-13 |
| side pairings are isometries; vertex cycles multiply to the identity | side residual 7.5e-28, cycles 2.0e-23 (160 bits) |
| Dirichlet polygon has area 272 pi and its sides pair up | difference -2.4e-14; endpoint mismatch 1.8e-27; 822 sides |
| mesh convergence: the 7140 distances between the 120 centers, coarse vs fine mesh | relative difference max 3.1e-2, median 8.2e-4 |
| the 5 real circles are closed geodesics (axis distance < 2e-2, length ratio within 1e-2, both shrinking on the fine mesh) | circles 1, 3, 4 pass (axis distance 1.5e-2, 7.2e-3, 9.8e-3); circle 2 fails (2.7e-2); circle 0 cannot be evaluated (see below) |
| charts reproduce the mesh values | chordal error of y: median 2.6e-6, max 2.2e-4, fraction above 1e-3: 0 |

## Known limitations

- Real circle 2: its largest distance from the geodesic is 2.7e-2, above the 2e-2 threshold, at the same place on both meshes
  (only about 3% of its points are above 1e-2). Not investigated further.
- Real circle 0 passes through the node y = 0. Near the node 1675 of its points could not be located on the mesh, and its
  holonomy has trace about 2e18, too large for the 200-bit axis computation; where the axis is undefined the viewer draws the
  circle as the geodesic through the ends of each run of points.
- Charts after the move onto the curve: 9808 of the 8 140 800 grid nodes moved more than 0.05 and were left where they were.
  Over all mesh vertices, with the best chart placement, the error of x is at most 2.1e-3 for 99.99% of vertices; the
  largest errors sit next to branch points and nodes.
- The mesh is not symmetric under complex conjugation, so the real-circle check is an independent test, but the picture is
  symmetric only up to the mesh error.
