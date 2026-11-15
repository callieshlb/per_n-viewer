#!/usr/bin/env python3
"""E92 explorer: export the data the browser needs into explorer/data.json."""
import os
import json
import sympy as sp
H = os.path.dirname(os.path.abspath(__file__)) + '/../'
x, y = sp.symbols('x y')
P = sp.sympify(open(os.path.dirname(os.path.abspath(__file__)) + '/../../inputs/E90-real-circles/per6_xy_poly.txt').read().replace('^', '**'))
dx = sp.Poly(P, x).degree()
coef = [[0]*8 for _ in range(dx+1)]           # coef[k][j]: coefficient of x^k y^j
for (kx, jy), c in sp.Poly(P, x, y).terms():
    coef[kx][jy] = int(c)
B = json.load(open(H + 'branch.json'))
cen = json.load(open(H + 'atlas_centers.json'))['centers']
G6 = json.load(open(os.path.dirname(os.path.abspath(__file__)) + '/../../inputs/E91-Gn-graph/G6.json'))
lab_of = {min(p): G6['labels'][i] for i, p in enumerate(G6['pairs'])}
pairs_of = {min(p): p for p in G6['pairs']}
centers = [dict(x=c['x'], y=c['y'], sheet=c['sheet'], label=lab_of[c['k']], angles=pairs_of[c['k']]) for c in cen]
out = dict(coef=coef, y0=B['y0'], x0=B['x0'], bps=B['bps'], perms=B['perms'], centers=centers)
json.dump(out, open(H + 'explorer/data.json', 'w'))
print('degree in x', dx, '| coefficient table', len(coef), 'x', len(coef[0]), '| centers', len(centers), '| branch points', len(B['bps']))

# punctures: census loop point at |g| = 1e-4 -> (x, y) = (a_5, a_3), sheet by continuation (as atlas_plot.py)
import re
import sys
import numpy as np
import mpmath as mp
sys.path.insert(0, H)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../../inputs/E90-real-circles')
import pern  # noqa: E402
from atlas_lib import label  # noqa: E402
census = json.load(open(os.path.dirname(os.path.abspath(__file__)) + '/../../inputs/E90-real-circles/data_n6_eps1e-4_dps60.json'))


def mpc_str(s):
    m = re.fullmatch(r'\(?\s*([-+]?[0-9.eE+-]+?)\s*([-+])\s*([0-9.eE+-]+)j\s*\)?', s.strip())
    return mp.mpc(mp.mpf(m.group(1)), mp.mpf(m.group(2)+m.group(3)))


mp.mp.dps = 40
pz = []
for zi, Pz in enumerate(census['punctures']):
    b, c = (mpc_str(s) for s in Pz['loop0'])
    a = pern.orbit(b, c, 6)
    ks = sorted(min(census['centers'][q['center']]['ks']) for q in Pz['polys'])
    pz.append(dict(x=complex(a[5]), y=complex(a[3]), name='/'.join(str(lab_of[k]) for k in ks),
                   real=bool(Pz['real_pts']), q=Pz['q']))
ty = np.array([p['y'] for p in pz])
ty = np.where(np.abs(ty-1) < 1e-13, 1+1e-13, ty)
X, _ = label(ty)
ch = lambda a, b: np.abs(a-b)/np.sqrt((1+np.abs(a)**2)*(1+np.abs(b)**2))
for p, row in zip(pz, X):
    p['sheet'] = int(ch(row, p['x']).argmin())
    p['x'] = [p['x'].real, p['x'].imag]
    p['y'] = [p['y'].real, p['y'].imag]
out['punctures'] = pz
json.dump(out, open(H + 'explorer/data.json', 'w'))
print('punctures', len(pz), [(p['name'], p['sheet']) for p in pz])
