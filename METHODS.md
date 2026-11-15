# Per_6(0): how the model was computed

## The curve

Per_6(0) is the curve of quadratic rational maps with a marked critical point of period 6. In the chart

    f(z) = (z - x)(z - w)/z^2,   w = (x + y - 1)/(x - 1),

the marked cycle is 0 -> inf -> 1 -> y -> f(y) -> x -> 0, and the curve is P(x, y) = 0 with deg_x P = 5.
The free critical point is omega = 2xw/(x + w). The compact curve has genus 6, 27 polynomial centers and 18 punctures.

## Steps

1. **Branched cover.** y : Per_6(0) -> sphere has degree 5. Branch points of y are computed, with the monodromy of a loop
   around each from a generic base point. The 5 sheets are labelled by continuing the roots x(y) along straight segments
   from the base point. (`atlas_lib.py`, `atlas_branch.py`)
2. **Mesh.** Each sheet is triangulated over a refined icosphere mesh of the y-sphere, with branch points as vertices, and
   the sheets are glued along the cuts by the monodromy. Result: one closed triangulated surface with 373 890 vertices,
   747 800 faces, Euler characteristic -10. (`surface3d/make_sphere_mesh.py`, `surface3d/edge_perms.js`, `surface3d/assemble.py`)
3. **Uniformization.** Edge lengths are chordal lengths on the y-sphere, so the mesh carries the complex structure of the
   curve. A discrete conformal change of metric (Bobenko-Pinkall-Springborn) is solved by damped Newton so that the angle
   sum is 2 pi at every vertex in the hyperbolic metric. (`surface3d/uniformize.py`)
4. **Development and fundamental polygon.** The surface, cut open to a disk, is developed into the Poincare disk. The side
   pairings generate the Fuchsian group; the Dirichlet polygon about a base point is the tile shown on the site
   (66 sides, 22 vertex cycles). (`surface3d/develop.py`, `surface3d/dirichlet_only.py`)
5. **Charts for the viewer.** (x, y) as functions on the disk are stored on 26 overlapping charts (160 x 160 grids), fitted
   to the mesh values by local holomorphic least squares. The viewer reads (x, y) at each pixel, polishes it on P = 0, and
   iterates the free critical point on the GPU. (`surface3d/chart_export.py`, `surface3d/conformal_export.py`)
6. **Overlays.** Polynomial centers, punctures, graph edges and the two real circles are located on the surface and
   exported. (`surface3d/overlay_export.py`, `surface3d/oval_check.py`, `explorer/export_data.py`, `surface3d/colour_export.py`)

## Checks

Each check was written down before it was run.

| check | result |
|---|---|
| Riemann-Hurwitz from the monodromy gives genus 6; product of the loops is the identity | ramification 20 at finite branch points, 0 at infinity; passes |
| the 27 centers lie on the curve, each in one sheet | distance to the labelled root <= 9e-14 |
| the 18 punctures are 18 distinct places | passes |
| angle sum 2 pi at every vertex | Newton residual 1.2e-13 (fine mesh) |
| total hyperbolic area = 20 pi (Gauss-Bonnet) | difference 2.8e-14 |
| side pairings are isometries matching the sides; vertex cycles multiply to the identity | side residual 4.8e-10, cycles 4.0e-8 (30-digit arithmetic) |
| Dirichlet polygon has area 20 pi and its sides pair up | difference -4.7e-10; 66 sides = 12g - 6 |
| mesh convergence: the 351 distances between the 27 centers, coarse vs fine mesh | relative difference: max 1.26e-2, median 3.37e-4 |
| the real circles are closed geodesics (they must be: complex conjugation is an anticonformal symmetry) | distance to the geodesic small and smaller on the fine mesh |

## Known limitations

- In double precision the side-pairing residual on the fine mesh is 1.5e-8, slightly above the 1e-8 that was set as the
  target; in 30-digit arithmetic it is 4.8e-10. The difference is round-off in the development.
- Neighbouring charts agree to better than 1e-3 on all but 5 of 6 806 overlapping pairs tested (worst 2.0e-2).
- The mesh is not symmetric under complex conjugation, so the real-circle check is an independent test, but also the picture
  is symmetric only up to the mesh error above.

## Files

See README.md.
