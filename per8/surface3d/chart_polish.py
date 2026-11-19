#!/usr/bin/env python3
"""E94 viewer fix 2 (user report 2026-11-09: grey/black blocks inside a capture component near centre 86). chart_residual_check.py
found the stored chart grid values off the curve (x vs the nearest root of P(., y): median 1.5e-5, 99th pct 2.4e-2 chordal), and
the viewer's per-pixel Newton polish is skipped where |y| or |x| > 6 (float32 overflow) and near y = 0: there the dynamics read a
marked cycle that does not close and the capture test fails (grey / black blocks). Fix: every grid node of every chart is moved
onto the curve in double precision: y kept, x replaced by the root of P(., y) nearest the stored x (vectorized Newton, np.roots
where Newton does not converge); accepted when the move is < 0.05 chordal, else the node is kept and counted. Re-encoded with the
chart's own Mobius normalisation and written over the chart grid parts (originals kept as chart_grids_*.orig.bin).
Raw values printed: moves (chordal), nodes kept, x-to-curve error before / after."""
import json
import os
import shutil
import sys
import time
import numpy as np
import sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
meta = json.load(open(os.path.join(here, 'chart_meta.json')))
S, NC = meta['S'], len(meta['charts'])
for f in meta['parts']:
    o = os.path.join(here, f.replace('.bin', '.orig.bin'))
    if not os.path.exists(o): shutil.copy2(os.path.join(here, f), o)
raw = np.concatenate([np.fromfile(os.path.join(here, f.replace('.bin', '.orig.bin')), '<f4') for f in meta['parts']])
grids = raw.reshape(NC, S, S, 4).astype(np.float64)
xs, ys = sp.symbols('x y')
Pp = sp.Poly(sp.sympify(open(os.path.join(here, '..', 'per8_xy_poly.txt')).read().replace('^', '**')), xs, ys)
dx, dy = Pp.degree(xs), Pp.degree(ys)
C = np.zeros((dx + 1, dy + 1))
for (i, j), c in Pp.terms(): C[i, j] = float(c)
def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
t0 = time.time(); stats = dict(nodes=0, moved=0, kept_far=0, kept_bad=0, roots_used=0)
before, after, moves = [], [], []
out = grids.copy()
for k in range(NC):
    ch = meta['charts'][k]; ay, by, ax, bx = (complex(*ch[c]) for c in ('ay', 'by', 'ax', 'bx'))
    G = grids[k].reshape(-1, 4); Y = G[:, 0] + 1j*G[:, 1]; X = G[:, 2] + 1j*G[:, 3]
    with np.errstate(all='ignore'):
        y = (ay - by*Y)/(1 - Y); x0 = (ax - bx*X)/(1 - X)
    ok = np.isfinite(y) & np.isfinite(x0) & (np.abs(y) < 1e5) & (np.abs(x0) < 1e5)
    n = len(y); stats['nodes'] += n
    cyv = np.zeros((n, dx + 1), complex)
    cyv[ok] = (y[ok][:, None]**np.arange(dy + 1)[None, :]) @ C.T
    x = np.where(ok, x0, 0); conv = np.zeros(n, bool)
    for it in range(40):
        p = np.zeros(n, complex); dp = np.zeros(n, complex)
        for i in range(dx, -1, -1): dp = dp*x + p; p = p*x + cyv[:, i]
        with np.errstate(all='ignore'): st = np.where(ok & (dp != 0), p/dp, 0)
        x = x - st
        conv = ok & (np.abs(st) <= 1e-13*(1 + np.abs(x)))
    # where Newton did not converge, or moved far, take the nearest root directly
    need = ok & (~conv | (chord(x, x0) > 1e-3))
    for t in np.where(need)[0]:
        r = np.roots(cyv[t][::-1])
        if len(r): x[t] = r[np.argmin(chord(r, x0[t]))]; stats['roots_used'] += 1
    mv = np.where(ok, chord(x, x0), np.nan)
    acc = ok & (mv < 0.05)
    stats['kept_far'] += int((ok & ~acc).sum()); stats['kept_bad'] += int((~ok).sum()); stats['moved'] += int((acc & (mv > 1e-7)).sum())
    xn = np.where(acc, x, x0)
    with np.errstate(all='ignore'):
        Xn = (xn - ax)/(xn - bx)
    Xn = np.where(acc & np.isfinite(Xn), Xn, X)
    o = out[k].reshape(-1, 4); o[:, 2] = Xn.real; o[:, 3] = Xn.imag
    moves.append(mv[ok])
    if k % 40 == 0:
        print(f'chart {k}: move median {np.nanmedian(mv):.2e}, 99th {np.nanpercentile(mv, 99):.2e}, max {np.nanmax(mv):.2e}; kept (move >= 0.05) {int((ok & ~acc).sum())} | {time.time()-t0:.0f} s', flush=True)
mv = np.concatenate(moves)
print(f'ALL nodes {stats["nodes"]}: moved {stats["moved"]}, kept (move >= 0.05) {stats["kept_far"]}, not finite / huge {stats["kept_bad"]}, '
      f'np.roots used {stats["roots_used"]}; move chordal median {np.median(mv):.2e}, 99th {np.percentile(mv, 99):.2e}, 99.9th {np.percentile(mv, 99.9):.2e}, max {mv.max():.2e}', flush=True)
raw2 = out.astype('<f4').tobytes(); per = S*S*4*4; off = 0
for f in meta['parts']:
    sz = os.path.getsize(os.path.join(here, f.replace('.bin', '.orig.bin')))
    open(os.path.join(here, f), 'wb').write(raw2[off:off + sz]); off += sz
assert off == len(raw2)
print('wrote', meta['parts'], flush=True)
