# Per_n(0) atlas

Interactive atlas of the curves Per_n(0) of quadratic rational maps, n = 1 ... 8 (Per_9 in preparation):
parameter space on the left (planes for n = 1..4, a flat torus for n = 5, uniformized hyperbolic models for n = 6, 7, 8),
the dynamical plane of the selected map on the right. Everything is computed live on the GPU (WebGL2).

**Viewing.** Once this repository is public, enable GitHub Pages (Settings → Pages → deploy from branch `main`, folder `/`)
and the site is at https://callieshlb.github.io/per_n-viewer/. Until then there is a private preview on claude.ai.
To run it locally instead: `python3 serve.py` (the pages fetch their data files, so opening index.html directly does not work).

Per_1 has external rays (parameter rays, and the dynamic rays of the same angles in the dynamical plane), with exact
landing points for rational angles.

Built from callieshlb/hurwitz, branch site-pern, runs/2026-11-09-E96-pern-site/make_site.py (see BUILD.txt).
