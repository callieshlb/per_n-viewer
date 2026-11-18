# Per_7(0): how the model was computed

## The curve

Per_7(0) is the curve of quadratic rational maps with a marked critical point of period 7. In the chart

    f(z) = (z - x)(z - w)/z^2,   w = (x + y - 1)/(x - 1),

the marked cycle is 0 -> inf -> 1 -> y -> a_4 -> a_5 -> x -> 0 (y = a_3, x = a_6), and the curve is P(x, y) = 0 with
deg_x P = deg_y P = 10. The free critical point is omega = 2xw/(x + w). The compact curve has genus 24, 63 polynomial
centers, 40 punctures and 4 real circles.

## Steps

1. **Branched cover.** y : Per_7(0) -> sphere has degree 10. Candidate branch points are the roots of the discriminant of P in x
   and of its leading coefficient (`disc.py`): 61 points. Over y = 0 the 10 roots form two clusters of 5 that are tangent to
   order 4, which double precision cannot separate for |y| < 0.3, so roots are continued in 256-bit arithmetic (`trackmp.py`).
   Monodromy of a loop around each candidate (`atlas_branch.py`): 54 simple branch points, two of cycle type (2,2,2,4) at
   y = 1/2 and y = 1, a node at y = 0 (trivial monodromy, two clusters of colliding roots), and 4 pole points.
   Riemann-Hurwitz: g = (66 - 20 + 2)/2 = 24.
2. **Centers and punctures on the sheets** (`atlas_points.py`): the 63 centers and the 40 puncture loop points of the census,
   located on the 10 sheets.
3. **Mesh of the y-sphere.** A refined icosphere with branch points as vertices, rings around each branch point, finer grids
   where branch points cluster, and no vertex within 0.03 of the node y = 0 (`surface3d/make_sphere_mesh.py`).
4. **Sheets and gluing.** Roots at every vertex and the sheet permutation along every edge, by 256-bit continuation
   (`surface3d/edge_perms.py`); the sheets glued into one closed surface (`surface3d/assemble.py`): 188 114 vertices,
   376 320 triangles, Euler characteristic -46.
5. **Uniformization.** Discrete conformal change of metric (Bobenko-Pinkall-Springborn) to a hyperbolic metric with angle sum
   2 pi at every vertex (`surface3d/uniformize.py`), then refined to 160-bit accuracy (`surface3d/refine.py`), because the
   development multiplies errors by up to e^23.
6. **Development and fundamental polygon.** Development into the Poincare disk in 160 bits (`surface3d/develop.py mp2`,
   `devmp.py`, `devmp2.py`); the side pairings generate the Fuchsian group; the Dirichlet polygon about a base point
   (282 = 12g - 6 sides, 94 vertex cycles) is the tile shown on the site (`surface3d/dirichlet_only.py`, `surface3d/fuchs.py`).
7. **Charts for the viewer.** (x, y) as functions on the disk on 112 overlapping charts (290 placements), fitted to the mesh by
   local holomorphic least squares (`surface3d/chart_export.py`); the chart map's edge pixels assigned to the nearest chart
   (`surface3d/chart_map_fill.py`). The viewer reads (x, y) per pixel, re-solves x on the curve, and iterates the free
   critical point on the GPU.
8. **Overlays.** Real circles (`surface3d/trace_ovals.py`, `surface3d/oval_check.py`), the graph G_7, centers and punctures
   (`surface3d/overlay_export.py`, `surface3d/conformal_export.py`, `surface3d/marker_info.py`).

## Checks

Each check was written down before it was run. A coarser mesh (not included) was computed alongside for the convergence checks.

| check | result |
|---|---|
| Riemann-Hurwitz from the monodromy; product of the loops = loop at infinity; monodromy group transitive | genus 24; both hold |
| mesh closed and orientable, Euler characteristic -46 (second route to the genus) | 0 closure failures, 0 orientation breaks; chi = -46 on both meshes |
| the 63 centers lie on the curve, each on one sheet, and are polynomial centers | max relative residual 7.9e-17; nearest root 4.6e-15, runner-up at least 2.1e14 times farther; all 63 attracting fixed points |
| the 40 punctures are 40 distinct places | passes |
| angle sum 2 pi at every vertex | 1.2e-13 in double precision; 8.5e-40 after refinement |
| total hyperbolic area = 4 pi (g - 1) = 92 pi | difference 2.9e-12 |
| side pairings are isometries; vertex cycles multiply to the identity | side residual 2.0e-29, cycles 1.6e-26 (160 bits) |
| Dirichlet polygon has area 92 pi and its sides pair up | difference 1.3e-11; endpoint mismatch 2.3e-11; 282 sides |
| mesh convergence: the 1953 distances between the 63 centers, coarse vs fine mesh | relative difference max 1.7e-2, median 5.6e-4 |
| the 4 real circles are closed geodesics, and closer to geodesics on the fine mesh than on the coarse one | distance to the axis 3.5e-3, 7.4e-3, 7.4e-3, 4.1e-3; length / geodesic length - 1: 2.0e-4, 9.8e-5, 6.3e-4, 2.9e-4. Circle 3 fails the second part (see below) |
| charts reproduce the mesh values | chordal error of y: median 6.5e-6, max 1.4e-4, fraction above 1e-3: 0 |

## Known limitations

- Real circle 3: its largest distance from the geodesic is larger on the fine mesh (4.1e-3) than on the coarse mesh
  (3.3e-3), so the registered "shrinks under refinement" test fails for it; its length ratio does improve, and all four
  circles meet the absolute thresholds (2e-2 and 1e-2). The cause has not been investigated.
- The double-precision root tracker in `explorer/per7_core.js` (used by the explorer and the viewer's clicks) cannot continue
  through the root clusters over y = 0: it fails for 20 of the 63 centers and 21 of the 40 punctures. The 256-bit tracker
  (`trackmp.py`) is used for everything that defines the model.
- The mesh is not symmetric under complex conjugation, so the real-circle check is an independent test, but the picture is
  symmetric only up to the mesh error.
