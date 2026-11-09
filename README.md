# Per_n(0) atlas

Interactive atlas of the curves Per_n(0) of quadratic rational maps, n = 1 ... 8 (Per_9 in preparation):
parameter space on the left (planes for n = 1..4, a flat torus for n = 5, uniformized hyperbolic models for n = 6, 7, 8),
the dynamical plane of the selected map on the right. Everything is computed live on the GPU (WebGL2).

**Private preview.** To look at it locally: `python3 serve.py` (opens http://127.0.0.1:8000/index.html).
The pages fetch their data files, so they need a server (opening index.html directly will not work).

Built from callieshlb/hurwitz, branch site-pern, runs/2026-11-09-E96-pern-site/make_site.py (see BUILD.txt).
