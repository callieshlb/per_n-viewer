#!/usr/bin/env python3
"""E95 vertex_x finishing step (replaces the tail of vertex_x.py, whose branch-vertex loop calls a 42-root solve for each of ~16 700
surface vertices over branch points and whose check (b) has no time cap): from the complete label table vx_labels.npz
- check (b): compare with the continuation roots of vx_checkb.npz (vertex_x_checkb.py), tolerance 1e-8 (1 + |x|);
- x per surface vertex: X[v, label] at normal vertices with |y| >= 0.1; NaN inside the masked disc and at the surface vertices
  over branch vertices (the chart fit uses their neighbours; no root solve at a multiple root).
Writes vertex_x.npz {x, y}. Exit 0 iff the BFS is complete and check (b) has 0 disagreements."""
import json, os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
os.chdir(here)
from meshloc import CS
M = json.load(open('sphere_mesh.json')); Y = np.array([complex(*v) for v in M['y']]); bp = np.array(M['bp'])
outer = set(int(v) for v in np.where(bp < 0)[0] if abs(Y[v]) >= 0.1)
z = np.load('vx_labels.npz', allow_pickle=True); seen = set(z['seen'].tolist()); front = list(z['front'])
X = np.load('vx_labels_polished.npz')['X']        # E95: roots Newton-polished at 512 bits (vertex_x_polish.py)
if os.environ.get('VX_FULL'):                     # E95: + the sheets over |y| < 0.1 (vertex_x_inner.py; 78 of 198 vertices by nearest-root fallback: approximate)
    zf = np.load('vx_labels_full.npz'); X = zf['X']; seen_full = set(zf['seen'].tolist())
    if os.environ.get('VX_FULL') == 'tracked': seen_full -= set(zf['approx'].tolist())   # only the exactly tracked inner vertices (the nearest-root ones failed check C)
print(f'labels: {len(seen)} of {len(outer)} outer normal vertices, BFS front {len(front)}, tree edges tracked {int(z["ntr"])}', flush=True)
complete = len(front) == 0 and seen == outer
cb = np.load('vx_checkb.npz'); bad = cmpd = 0; worst = 0.0
for v, R, st in zip(cb['v'], cb['R'], cb['st']):
    if st != 1 or int(v) not in seen: continue
    h = CS.h_at(complex(Y[v])); cmpd += 1
    dv = max(abs(X[v, h[j]] - R[j])/(1 + abs(R[j])) for j in range(42)); worst = max(worst, dv)
    if dv > 1e-8: bad += 1
print(f'check (b): {cmpd} vertices compared with the continuation from y0 (of 200: failed {int((cb["st"] == 2).sum())}, capped {int((cb["st"] == 3).sum())}, '
      f'not yet run {int((cb["st"] == 0).sum())}): label disagreements {bad}; max relative difference {worst:.2e}', flush=True)
S = np.array(json.load(open('surface_mesh.json'))['V'])
v, l = S[:, 0], S[:, 1]
xs = np.full(len(S), np.nan + 0j); ys = Y[v]
okm = (l >= 0) & np.isin(v, np.array(sorted(seen_full if os.environ.get('VX_FULL') else seen)))
xs[okm] = X[v[okm], l[okm]]
print(f'surface vertices {len(S)}: x finite {int(np.isfinite(xs).sum())}; NaN {int((~np.isfinite(xs)).sum())} = over branch vertices {int((l < 0).sum())} + inside |y| < 0.1 or unlabelled '
      f'{int(((l >= 0) & ~okm).sum())}; max |x| {np.nanmax(np.abs(xs)):.3e}', flush=True)
np.savez('vertex_x.npz', x=xs, y=ys)
sys.exit(0 if complete and bad == 0 and cmpd > 0 else 1)
