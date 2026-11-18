# Per_6(0): how the model was computed

Version 2 (2026-11-18): the mesh is refined at the punctures, and the colouring near the puncture 6/58 uses an extra chart.
Version 1 (the first, coarser mesh near the punctures) is in the history of this branch.

## The curve

Per_6(0) is the curve of quadratic rational maps with a marked critical point of period 6. In the chart

    f(z) = (z - x)(z - w)/z^2,   w = (x + y - 1)/(x - 1),

the marked cycle is 0 -> inf -> 1 -> y -> f(y) -> x -> 0, and the curve is P(x, y) = 0 with deg_x P = 5.
The free critical point is omega = 2xw/(x + w). The compact curve has genus 6, 27 polynomial centers and 18 punctures.

## Steps

1. **Branched cover.** y : Per_6(0) -> sphere has degree 5. Branch points of y with the monodromy of a loop around each, from a
   generic base point; sheets labelled by continuing the roots x(y) along straight segments from the base point.
   (`atlas_lib.py`, `atlas_branch.py`)
2. **Mesh of the y-sphere.** A refined icosphere with branch points as vertices, plus geometric rings around the seven values
   of y over which the punctures lie (1/2, 1, 0, 2/3, (sqrt5 - 1)/2, -(sqrt5 + 1)/2, infinity), shrinking to radius 1e-9 at the
   cone points 1/2 and 1, and to 1e-4 to 5e-4 at the others. (`surface3d/make_sphere_mesh.py`)
3. **Sheets and gluing.** The 5 roots at every vertex and the sheet permutation along every edge (`surface3d/edge_perms.js`),
   recomputed in 40-digit arithmetic near y = 0, 1/2, 1 where roots cluster (`surface3d/edge_perms_mp.py`); the sheets glued
   into one closed surface (`surface3d/assemble.py`): 472 800 vertices, 945 620 triangles, Euler characteristic -10.
4. **Uniformization.** Edge lengths are chordal lengths on the y-sphere, so the mesh carries the complex structure of the curve.
   A discrete conformal change of metric (Bobenko-Pinkall-Springborn) makes the angle sum 2 pi at every vertex in the
   hyperbolic metric. (`surface3d/uniformize.py`)
5. **Development and fundamental polygon.** The cut-open surface is developed into the Poincare disk in 30-digit arithmetic;
   the side pairings generate the Fuchsian group; the Dirichlet polygon about a base point (66 sides, 22 vertex cycles) is
   the tile shown on the site. (`surface3d/develop.py`, `surface3d/devmp.py`, `surface3d/dirichlet_only.py`)
6. **Charts for the viewer.** (x, y) as functions on the disk on 24 overlapping charts (160 x 160 grids), fitted to the mesh
   by local holomorphic least squares (`surface3d/chart_export.py`); pixels of the chart map at the edge of the disk assigned to
   the nearest chart (`surface3d/chart_map_fill.py`). The viewer reads (x, y) per pixel, polishes x on P = 0, and iterates
   the free critical point on the GPU.
7. **Cone point 6/58.** The 6/58 puncture is a cone point of y. Within hyperbolic distance 0.15 of it, a 25th chart stores
   the offsets y - y_c, x - x_c from the cone values, computed from float64 grids, so the single-precision GPU arithmetic does
   not cancel there. (`cone/chart_export64.py`, `cone/cone_data.py`)
8. **Overlays.** Polynomial centers, punctures, the graph G_6 and the two real circles located on the surface and exported.
   (`surface3d/oval_check.py`, `surface3d/overlay_export.py`, `surface3d/conformal_export.py`)

## Checks

Each check was written down before it was run.

| check | result |
|---|---|
| mesh closed and orientable, Euler characteristic -10 | 0 closure failures, 0 orientation breaks; chi = -10 |
| sheet permutations near y = 0, 1/2, 1 (40 digits) | 18 648 edges, 0 failures, 0 non-permutation maps |
| angle sum 2 pi at every vertex | max error 2.1e-13 |
| total hyperbolic area = 20 pi | difference 5.7e-13 |
| side pairings are isometries; vertex cycles multiply to the identity | side residual 3.0e-9, cycles 5.3e-8 (30 digits) |
| Dirichlet polygon has area 20 pi and its sides pair up | difference 9.0e-10; endpoint mismatch 3.8e-9; 66 sides = 12g - 6 |
| the 351 distances between the 27 centers agree with the earlier (version 1) mesh | relative difference max 2.2e-3, median 6.4e-5 |
| the real circles are closed geodesics (complex conjugation is an anticonformal symmetry) | distance to the axis 1.6e-3 and 1.7e-3; length / geodesic length - 1: 1.0e-4 and 1.1e-4 |
| near each puncture (hyperbolic 0.1) all triangles have diameter <= 0.025 | 17 of 18 punctures; at puncture 30 the largest is 0.0255 |
| every puncture has a triangle centre within 0.02 | 18 of 18 (nearest at most 1.1e-3) |

## Known limitations

- At puncture 30 one triangle is slightly larger than the target (0.0255 against 0.025).
- At puncture 32 (also a cone point) a small zone (scattered points within about 0.012 in the hyperbolic metric) is still not resolved:
  the single-precision chart storage cannot separate the values there. Puncture 28 (a node) shows a noisy tuft.
- Next to 6/58, at the viewer's default 600 iterations, a black band is mostly "undecided" (slow, parabolic convergence of the
  free critical point), not a boundary; at 60 000 iterations 1-3% of the pixels near the cone stay undecided.
- At puncture 2/62 there are small jagged patches along some overlay lines; not investigated.
- The mesh is not symmetric under complex conjugation, so the real-circle check is an independent test, but the picture is
  symmetric only up to the mesh error above.
