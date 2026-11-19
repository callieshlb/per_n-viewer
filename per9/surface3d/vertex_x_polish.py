#!/usr/bin/env python3
"""E95: Newton polish (512 bits, trackmp.Curve.newton) of every labelled root in vx_labels.npz: the vertex roots were the midpoints
of flint's isolating balls (radius up to 1e-2; values off by up to ~5e-6), which check (b) exposed (labels identical, values not).
Each polished root must stay within its own ball: moved by less than half the distance to the nearest other root at the vertex
(counted otherwise). Writes vx_labels_polished.npz {X, seen}. 4 workers, checkpoint per block."""
import json, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
os.chdir(here)
from meshloc import CURVE as C
from gmpy2 import mpc
M = json.load(open('sphere_mesh.json')); Y = np.array([complex(*v) for v in M['y']])
z = np.load('vx_labels.npz', allow_pickle=True); X = z['X']; seen = np.array(sorted(z['seen'].tolist()))
assert len(z['front']) == 0
BLK = 2000; os.makedirs('vx_polish', exist_ok=True)
def blk(b):
    fn = f'vx_polish/b{b:04d}.npz'
    if os.path.exists(fn): return 0, 0.0
    vs = seen[b*BLK:(b + 1)*BLK]; out = np.empty((len(vs), 42), complex); far = 0; mv = 0.0
    for i, v in enumerate(vs):
        cc = C.cy(mpc(complex(Y[v])))[0]; x0 = X[v]
        xp = np.array([complex(C.newton(cc, mpc(complex(x)))[0]) for x in x0])
        sep = np.abs(x0[:, None] - x0[None, :]) + np.eye(42)*1e300
        d = np.abs(xp - x0); far += int((d > 0.5*sep.min(1)).sum()); mv = max(mv, float((d/(1 + np.abs(x0))).max()))
        out[i] = xp
    np.savez(fn + '.tmp.npz', X=out, far=far, mv=mv); os.replace(fn + '.tmp.npz', fn); return far, mv
if __name__ == '__main__':
    from multiprocessing import Pool
    t0 = time.time(); nb = (len(seen) + BLK - 1)//BLK
    with Pool(int(os.environ.get('NW', '4'))) as pool: pool.map(blk, range(nb), chunksize=1)
    far = 0; mv = 0.0; XP = X.copy()
    for b in range(nb):
        zz = np.load(f'vx_polish/b{b:04d}.npz'); XP[seen[b*BLK:(b + 1)*BLK]] = zz['X']; far += int(zz['far']); mv = max(mv, float(zz['mv']))
    print(f'polished {len(seen)} x 42 roots: max relative move {mv:.2e}; roots that moved more than half the distance to their nearest neighbour {far}; '
          f'finite {bool(np.isfinite(XP[seen]).all())}; {time.time()-t0:.0f} s', flush=True)
    np.savez('vx_labels_polished.npz', X=XP, seen=seen)
    sys.exit(0 if far == 0 else 1)
