#!/usr/bin/env python3
"""E93 A2/A3 and explorer data (generalizes E92 atlas_render.py's A2, atlas_plot.py's A3 and explorer/export_data.py).
Centres: census c (z^2 + c, period-7 critical point) -> pern.center_to_chart -> pern.orbit -> (x, y) = (a_6, a_3).
Punctures: census loop point loop0 -> pern.orbit -> (a_6, a_3) (as E92: the loop point near the puncture, |g| = eps).
Atlas sheet of (y, x): continue the 10 base roots X0 (branch.json, lexsorted at y0) from y0 to y along the straight segment
(node detour radius 0.025 round y = 0), 256-bit (trackmp.py); sheet = index of the continued root nearest x.
A2: |P(x, y)| / sum |terms| < 1e-10; distance to the nearest root < 1e-8 with runner-up > 1e3 times farther; fate of the free
critical point = attracting cycle of period 1 (checked by explorer/test_core.js with the explorer core).
A3: 40 distinct places (sheet, y rounded to 1e-2).
Writes atlas_points.json and explorer/data.json."""
import json
import os
import re
import sys
import time
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
E90 = os.path.join(HERE, '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import pern  # noqa: E402
import trackmp as tm  # noqa: E402
from gmpy2 import mpc  # noqa: E402
B = json.load(open(os.path.join(HERE, 'branch.json')))
BPS = [complex(*b) for b in B['bps']]
NODES = [(complex(*b), 0.025) for b, k in zip(B['bps'], B['kinds']) if k == 'node']
C = tm.Curve(specials=BPS); n_sh = C.n
y0 = complex(*B['y0'])
X0 = [mpc(complex(*v)) for v in B['x0']]
c0 = C.cy(mpc(y0))[0]
X0 = [C.newton(c0, x)[0] for x in X0]
census = json.load(open(os.path.join(E90, 'data_n7_eps1e-3_dps60.json')))
G7 = json.load(open(os.path.join(HERE, '..', 'inputs', 'E91-Gn-graph', 'G7.json')))
lab_of = {min(p): G7['labels'][i] for i, p in enumerate(G7['pairs'])}
pairs_of = {min(p): p for p in G7['pairs']}
def mpc_str(s):
    m = re.fullmatch(r'\(?\s*([-+]?[0-9.eE+-]+?)\s*([-+])\s*([0-9.eE+-]+)j\s*\)?', s.strip())
    return mp.mpc(mp.mpf(m.group(1)), mp.mpf(m.group(2)+m.group(3)))
def atlas(y, x):
    st = {}
    xe = tm.track_line(C, X0, y0, y, NODES, 0.15, st)
    if xe is None:
        return None, None, None, st
    d = sorted((float(abs(xx - x)), k) for k, xx in enumerate(xe))
    return d[0][1], d[0][0], d[1][0], st
def onP(y, x):
    """|P| / sum |terms| in 256-bit arithmetic."""
    c = C.cy(mpc(y))[0]; X = mpc(x)
    val = sum(ck*X**k for k, ck in enumerate(c))
    tot = sum(abs(a)*abs(mpc(y))**j*abs(X)**k for k, row in enumerate(C.C) for j, a in enumerate(row))
    return float(abs(val)/tot)
t0 = time.time()
mp.mp.dps = 40
cen = []
for e in census['centers']:
    b, c, _ = pern.center_to_chart(mp.mpc(mp.mpf(e['c'][0]), mp.mpf(e['c'][1])))
    a = pern.orbit(b, c, 7)
    x, y = complex(a[6]), complex(a[3])
    k = min(e['ks'])
    rts = sorted(float(abs(r - mpc(x))) for r in C.roots(y))
    s, d1, d2, st = atlas(y, x)
    cen.append(dict(x=[x.real, x.imag], y=[y.real, y.imag], k=k, label=lab_of[k], angles=pairs_of[k], sheet=s,
                    relP=onP(y, x), droot=rts[0], ratio=rts[1]/max(rts[0], 1e-300), atlas_d=d1, atlas_d2=d2, real=bool(e['real'])))
print(f'A2 centres: {len(cen)} | max |P|/sum|terms| {max(c["relP"] for c in cen):.2e} | max distance to nearest root '
      f'{max(c["droot"] for c in cen):.2e} | min runner-up/nearest ratio {min(c["ratio"] for c in cen):.2e} | continuation failures '
      f'{sum(c["sheet"] is None for c in cen)} | max distance to the continued root {max(c["atlas_d"] or 9 for c in cen):.2e}, '
      f'min runner-up {min(c["atlas_d2"] or 9 for c in cen):.2e} | centres per atlas sheet '
      f'{[sum(1 for c in cen if c["sheet"] == s) for s in range(n_sh)]} | {time.time()-t0:.1f} s', flush=True)
mp.mp.dps = 40
pz = []
for zi, Pz in enumerate(census['punctures']):
    b, c = (mpc_str(s) for s in Pz['loop0'])
    a = pern.orbit(b, c, 7)
    x, y = complex(a[6]), complex(a[3])
    ks = sorted(min(census['centers'][q['center']]['ks']) for q in Pz['polys'])
    rts = sorted(float(abs(r - mpc(x))) for r in C.roots(y))
    s, d1, d2, st = atlas(y, x)
    pz.append(dict(x=[x.real, x.imag], y=[y.real, y.imag], name='/'.join(str(lab_of[k]) for k in ks), real=bool(Pz['real_pts']),
                   q=Pz['q'], sheet=s, relP=onP(y, x), droot=rts[0], ratio=rts[1]/max(rts[0], 1e-300), atlas_d=d1, atlas_d2=d2))
places = {(p['sheet'], round(p['y'][0], 2), round(p['y'][1], 2)) for p in pz}
print(f'A3 punctures: {len(pz)} loop points | continuation failures {sum(p["sheet"] is None for p in pz)} | distinct places (sheet, y rounded '
      f'1e-2): {len(places)} | max |P|/sum|terms| {max(p["relP"] for p in pz):.2e} | max distance to nearest root '
      f'{max(p["droot"] for p in pz):.2e} | min runner-up/nearest {min(p["ratio"] for p in pz):.2e}', flush=True)
for p in pz:
    print(f'   puncture {p["name"]:>9s} q {p["q"]} y {complex(*p["y"]):.6f} x {complex(*p["x"]):.6f} sheet {p["sheet"]} '
          f'd {p["atlas_d"]:.1e} runner-up {p["atlas_d2"]:.1e}')
json.dump(dict(centers=cen, punctures=pz), open(os.path.join(HERE, 'atlas_points.json'), 'w'), indent=1)
# explorer data (format of E92 explorer/data.json + n, leadQ, leadPow, nodes, kinds)
coef = [[int(c) for c in row] for row in C.coef_int]
out = dict(n=7, coef=coef, leadQ=[1, -2, -2, 3, 1], leadPow=2, y0=B['y0'], x0=[[complex(x).real, complex(x).imag] for x in X0],
           bps=B['bps'], kinds=B['kinds'], perms=B['perms'], nodes=[[0.0, 0.0, 0.2]],
           centers=[dict(x=c['x'], y=c['y'], sheet=c['sheet'], label=c['label'], angles=c['angles']) for c in cen],
           punctures=[dict(x=p['x'], y=p['y'], sheet=p['sheet'], name=p['name'], real=p['real'], q=p['q']) for p in pz])
os.makedirs(os.path.join(HERE, 'explorer'), exist_ok=True)
json.dump(out, open(os.path.join(HERE, 'explorer', 'data.json'), 'w'))
print('wrote atlas_points.json, explorer/data.json')
