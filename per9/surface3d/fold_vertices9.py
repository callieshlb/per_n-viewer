#!/usr/bin/env python3
"""E95 viewer step: every surface vertex of the mesh folded into the certified Dirichlet polygon D (develop.npz), as E94's
chart_export.py does (P = fold(z[first])), with an active-set fold (4.07 M points x 2250 faces). Writes P_folded.npy (complex128
per surface vertex) and prints raw checks. usage: fold_vertices9.py [mesh_dir]"""
import os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
d = sys.argv[1] if len(sys.argv) > 1 else here
t0 = time.time()
D = np.load(os.path.join(d, 'develop.npz')); z, orig, FG = D['z'], D['orig'], D['face_g']
U, T, _ = fu.halfplanes(FG); Gi = fu.su11(np.linalg.inv(FG))
Ur = np.stack([U.real, U.imag])                    # 2 x F
nV = int(orig.max()) + 1
first = np.full(nV, -1); first[orig[::-1]] = np.arange(len(orig))[::-1]
assert (first >= 0).all()
w = z[first].copy()
print(f'surface vertices {nV}; developed radius of the chosen lifts: max {fu.hd0(w).max():.3f}', flush=True)
active = np.arange(nV); CH = 40000; it = 0
while len(active) and it < 400:
    nxt = []
    for s0 in range(0, len(active), CH):
        idx = active[s0:s0 + CH]; k = fu.to_klein(w[idx])
        v = np.stack([k.real, k.imag], 1) @ Ur - T[None, :]
        a = np.argmax(v, 1); out = v[np.arange(len(idx)), a] > 1e-13
        if out.any():
            g = Gi[a[out]]; ww = w[idx[out]]
            w[idx[out]] = (g[:, 0, 0]*ww + g[:, 0, 1])/(g[:, 1, 0]*ww + g[:, 1, 1])
            nxt.append(idx[out])
    active = np.concatenate(nxt) if nxt else np.array([], int); it += 1
    print(f'   fold pass {it}: still outside {len(active)}, {time.time()-t0:.0f} s', flush=True)
k = fu.to_klein(w); worst = -np.inf
for s0 in range(0, nV, CH):
    kk = k[s0:s0 + CH]; worst = max(worst, float((np.stack([kk.real, kk.imag], 1) @ Ur - T[None, :]).max()))
print(f'folded {nV} vertices in {it} passes; not converged {len(active)}; max half-plane violation {worst:.2e}; finite {bool(np.isfinite(w).all())}; '
      f'hyperbolic radius in D: median {np.median(fu.hd0(w)):.3f}, 90% {np.percentile(fu.hd0(w), 90):.3f}, 99% {np.percentile(fu.hd0(w), 99):.3f}, max {fu.hd0(w).max():.3f}; '
      f'fraction beyond 5.5: {(fu.hd0(w) > 5.5).mean():.4f}, beyond 6.9: {(fu.hd0(w) > 6.9).mean():.4f}; {time.time()-t0:.0f} s', flush=True)
np.save(os.path.join(d, 'P_folded.npy'), w)
sys.exit(0 if len(active) == 0 and worst < 1e-9 else 1)
