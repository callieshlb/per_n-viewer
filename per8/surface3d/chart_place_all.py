#!/usr/bin/env python3
"""E94 viewer fix (user report 2026-11-09: straight-edged glitches). chart_dense_check.py found 83 551 of 844 924 vertices read
through no valid chart: the 1024^2 base-frame chart map leaves texels unassigned when the texel centre is outside D (a band
along every edge of D), and near the rim one texel spans ~0.6 hyperbolic, more than a chart can cover.
Fix, part 1 (this script): export ALL chart placements covering D, i.e. every image h(c_k) of a chart centre under
H = {I} u F u F.F (as chart_export.py) within rD + R_COV + 0.3 of 0, as chart_meta.json['cand'] = [[k, M (8 floats), p (2)], ...]
with M = Mq(c_k)^-1 h^-1 (base frame -> chart coordinate u) and p = h(c_k). Part 2 (make_live.py): the viewer picks, per view,
the placements near the view centre (tiles found by BFS, double precision) and each pixel uses the one it is deepest in.
Check printed: for every surface vertex (folded into D), the best placement in 'cand' (smallest |u|) and its chart value
error in y and x (as chart_dense_check.py); coverage = max over vertices of min |u| (must be < 0.97 RU)."""
import base64
import json
import os
import sys
import numpy as np
from scipy.spatial import cKDTree
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from fuchs import Mq, mob_ap, hd0
from meshloc import load_D
import chart_dense_check as cdc          # reuses its loaded data and read(); its own check runs on import (prints)

meta = json.load(open('chart_meta.json'))
res, D = load_D('.')
rD = float(hd0(fu.to_poinc(res['verts'])).max()); R_COV = 1.2
Fmp = fu.load_face_mp(D)
Fall = fu.mp_dedupe(Fmp + [fu.mp_inv(e) for e in Fmp])
Fd = np.array([fu.mp_to_np(e) for e in Fall])
cent = np.array([complex(*c['c']) for c in meta['charts']])
# images of the chart centres under H, screened in double, elements formed in 200 bits only for kept images
cand = {}
IDC = fu.gmpy2.cosh(fu.mpfr('5e-4'))
def keep(k, h):
    hn = fu.mp_to_np(h); p = mob_ap(hn, cent[k])
    if hd0(np.array([p]))[0] > rD + R_COV + 0.3: return
    key = (k, round(p.real*1e7), round(p.imag*1e7))
    if key in cand: return
    Mp = np.linalg.inv(Mq(cent[k])) @ np.linalg.inv(hn); Mp = Mp/np.sqrt(np.linalg.det(Mp))
    cand[key] = [int(k)] + [float(v) for z_ in (Mp[0, 0], Mp[0, 1], Mp[1, 0], Mp[1, 1]) for v in (z_.real, z_.imag)] + [p.real, p.imag]
I2 = (fu.mpc(1), fu.mpc(0))
for k in range(len(cent)):
    keep(k, I2)
    for e in Fall: keep(k, e)
for s0 in range(0, len(Fall), 50):
    P2 = fu.su11(np.einsum('nij,mjk->nmik', Fd[s0:s0 + 50], Fd).reshape(-1, 2, 2))
    for k in range(len(cent)):
        im = (P2[:, 0, 0]*cent[k] + P2[:, 0, 1])/(P2[:, 1, 0]*cent[k] + P2[:, 1, 1])
        for q in np.where(hd0(im) <= rD + R_COV + 0.35)[0]:
            a_, b_ = divmod(int(q), len(Fall))
            h = fu.mp_mul(Fall[s0 + a_], Fall[b_])
            if abs(h[0]) < IDC: continue
            keep(k, h)
C = list(cand.values())
print(f'cand placements: {len(C)} (charts {len(cent)}, rD {rD:.4f})', flush=True)
# check: every vertex has a placement with |u| < 0.97 RU, and its chart value is right
CP = np.array([complex(c[9], c[10]) for c in C]); CK = np.array([c[0] for c in C])
CM = np.array([[[c[1] + 1j*c[2], c[3] + 1j*c[4]], [c[5] + 1j*c[6], c[7] + 1j*c[8]]] for c in C])
X3 = np.stack([CP.real, CP.imag], 1)
P = cdc.P; tree = cKDTree(X3)
_, nb = tree.query(np.stack([P.real, P.imag], 1), 24)
best_u = np.full(len(P), np.inf); best_j = np.zeros(len(P), int)
for t in range(nb.shape[1]):
    j = nb[:, t]; M = CM[j]; u = (M[:, 0, 0]*P + M[:, 0, 1])/(M[:, 1, 0]*P + M[:, 1, 1]); a = np.abs(u)
    b = a < best_u; best_u[b] = a[b]; best_j[b] = j[b]
miss = best_u > 0.97*cdc.RU
print(f'vertices whose best of the 24 Euclidean-nearest placements has |u| > 0.97 RU: {int(miss.sum())} of {len(P)}; max min|u| / RU {best_u.max()/cdc.RU:.4f}', flush=True)
ey = np.zeros(len(P)); ex = np.zeros(len(P))
M = CM[best_j]; u = (M[:, 0, 0]*P + M[:, 0, 1])/(M[:, 1, 0]*P + M[:, 1, 1]); kk = CK[best_j]
for k in np.unique(kk):
    m = kk == k; y_, x_ = cdc.read(k, u[m]); ey[m] = cdc.chord(y_, cdc.yv[m]); ex[m] = cdc.chord(x_, cdc.xv[m])
for name, e in (('y', ey), ('x', ex)):
    print(f'  best-placement chordal {name} error: median {np.median(e):.2e}, 99th pct {np.percentile(e, 99):.2e}, 99.99th {np.percentile(e, 99.99):.2e}, '
          f'max {e.max():.2e}, count > 1e-3 {int((e > 1e-3).sum())}', flush=True)
meta['cand'] = C
json.dump(meta, open('chart_meta.json', 'w'))
print('wrote chart_meta.json with cand', flush=True)
