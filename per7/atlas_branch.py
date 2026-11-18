#!/usr/bin/env python3
"""E93 step A (generalizes E92 atlas_branch.py to n = 7): branch points of y over P (per7_xy_poly.txt), lasso monodromy at a
base point y0, monodromy at infinity, Riemann-Hurwitz genus (A1). Root continuation in 256-bit arithmetic (trackmp.py);
run 1 with E92's np.roots matching FAILED near y = 0 (log_branch_run1_nproots_FAILED.txt: order-4 tangency there).
Each lasso is computed twice (step factor 0.15 and 0.05; must agree). Kinds: 'branch' (non-trivial lasso), 'node'
(trivial lasso, roots collide as y -> b), 'pole' (trivial lasso, one root -> x = inf, no collision).
Writes branch.json. Raw values only."""
import cmath
import json
import math
import os
import sys
import time
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trackmp import Curve, match, tofloat
x, y = sp.symbols('x y')
P = sp.sympify(open(os.path.join(HERE, '..', 'inputs', 'E90-real-circles', 'per7_xy_poly.txt')).read().replace('^', '**'))
lead = sp.Poly(P, x).all_coeffs()[0]
cand = []
for line in open(os.path.join(HERE, 'disc_factors.txt')):
    e, f = line.split(' ', 1); f = sp.sympify(f)
    for r in sp.Poly(f, y).nroots(n=40, maxsteps=500):
        cand.append((complex(r), f'disc factor deg {sp.Poly(f, y).degree()} mult {e}'))
for f, e in sp.factor_list(lead)[1]:
    for r in sp.Poly(f, y).nroots(n=40, maxsteps=500):
        r = complex(r)
        if all(abs(r - o) > 1e-10 for o, _ in cand):
            cand.append((r, f'leading coefficient factor deg {sp.Poly(f, y).degree()} only'))
bps = [c for c, _ in cand]; src = [s for _, s in cand]
C = Curve(specials=bps); d = C.n
dmin = min(abs(a - b) for i, a in enumerate(bps) for b in bps[i+1:])
print('degree in x', d, '| candidates', len(bps), '| min pairwise distance', f'{dmin:.3e}', '| max modulus', f'{max(abs(b) for b in bps):.4f}')
y0 = complex(0.31, -0.47)
print('base point y0', y0, '| distance to nearest candidate', f'{min(abs(y0 - b) for b in bps):.4f}', flush=True)
X0 = C.roots(y0)
f0 = np.array(tofloat(X0)); o = np.lexsort((f0.imag, f0.real)); X0 = [X0[i] for i in o]
def perm_of(xe):
    p, w = match(xe, X0)
    assert sorted(p) == list(range(d)), ('not a permutation', p)
    return p, w
def lasso(b, rad, hfac, m=96):
    st = b + rad*(y0-b)/abs(y0-b)
    th0 = cmath.phase(st-b)
    circ = [b+rad*cmath.exp(1j*(th0+2*math.pi*t/m)) for t in range(m+1)]
    stats = {}
    xe = C.track_path(X0, [y0] + circ + [y0], hfac, stats)
    if xe is None: return None, None, stats
    p, w = perm_of(xe); return p, w, stats
t0 = time.time(); perms, rads, worst, disagree = [], [], 0.0, 0
for k, b in enumerate(bps):
    rad = min(min(abs(b - o) for o in bps if o != b)/3, 0.2)
    p1, w1, s1 = lasso(b, rad, 0.15); p2, w2, s2 = lasso(b, rad, 0.05)
    if p1 is None or p2 is None:
        print('TRACKING FAILED at candidate', k, b, s1, s2); sys.exit(1)
    if p1 != p2: disagree += 1; print('lasso disagreement at', b, p1, p2)
    worst = max(worst, w1, w2); perms.append(p1); rads.append(rad)
print(f'lassos done in {time.time()-t0:.1f} s; step-factor 0.15 vs 0.05 disagreements {disagree}; worst match ratio {worst:.2e}', flush=True)
R = 4*max(abs(b) for b in bps) + 4
circ = [y0 + R*cmath.exp(2j*math.pi*t/4000) for t in range(4001)]
pinf_list = []
for hfac in (0.15, 0.05):
    xe = C.track_path(X0, [y0] + circ + [y0], hfac)
    pinf_list.append(perm_of(xe)[0])
pinf = pinf_list[0]
print('loop at infinity (radius', R, ') step factors 0.15 / 0.05 agree:', pinf_list[0] == pinf_list[1])
def cycles(p):
    seen, cs = set(), []
    for i in range(d):
        if i not in seen:
            c = []; j = i
            while j not in seen:
                seen.add(j); c.append(j); j = p[j]
            cs.append(c)
    return cs
Rfin = sum(d-len(cycles(p)) for p in perms)
Rinf = d-len(cycles(pinf))
g = (Rfin+Rinf-2*d+2)/2
order = sorted(range(len(bps)), key=lambda k: cmath.phase(bps[k]-y0) % (2*math.pi))
comp = list(range(d))
for k in order:
    comp = [perms[k][comp[i]] for i in range(d)]
reach = {0}; fr = [0]
while fr:
    i = fr.pop()
    for p in perms:
        if p[i] not in reach: reach.add(p[i]); fr.append(p[i])
def collide(b):
    yy = b + 1e-7*(y0-b)/abs(y0-b)
    r = np.array(tofloat(C.roots(yy)))
    dd = np.abs(r[:, None] - r[None, :])/np.sqrt((1+np.abs(r[:, None])**2)*(1+np.abs(r[None, :])**2))
    return int((dd[np.triu_indices(d, 1)] < 1e-3).sum()), float(np.abs(r).max())
kinds, coll = [], []
for b, p in zip(bps, perms):
    npair, xmax = collide(b); coll.append((npair, xmax))
    kinds.append('branch' if p != list(range(d)) else ('node' if npair > 0 else 'pole'))
print('non-trivial lassos', kinds.count('branch'), '| nodes (trivial lasso, colliding roots)', kinds.count('node'),
      '| pole points (trivial lasso, no collision)', kinds.count('pole'))
for b, p, s, k, cl in zip(bps, perms, src, kinds, coll):
    cyc = [c for c in cycles(p) if len(c) > 1]
    if k != 'branch' or len(cyc) != 1 or len(cyc[0]) != 2:
        print('   ', k, f'y = {b.real:.12f}{b.imag:+.12f}i', '|', s, '| cycles (len>1)', cyc, '| colliding pairs at 1e-7, max|x|', cl)
print('cycle-type histogram of non-trivial lassos:',
      {str(t): sum(1 for p in perms if tuple(sorted(len(c) for c in cycles(p) if len(c) > 1)) == t and p != list(range(d)))
       for t in sorted({tuple(sorted(len(c) for c in cycles(p) if len(c) > 1)) for p in perms if p != list(range(d))})})
print('ramification finite', Rfin, '| at infinity', Rinf, '| loop at infinity perm', pinf, '| cycles', [c for c in cycles(pinf) if len(c) > 1])
print('product of lassos (counterclockwise order about y0)', comp, '| equals loop at infinity:', comp == pinf)
print('monodromy group transitive:', len(reach) == d)
print('Riemann-Hurwitz genus', g)
json.dump(dict(y0=[y0.real, y0.imag], x0=[[complex(v).real, complex(v).imag] for v in X0], R=R, bps=[[b.real, b.imag] for b in bps],
               perms=perms, kinds=kinds, src=src, pinf=pinf, genus=g, Rfin=Rfin, Rinf=Rinf, lasso_radius=rads),
          open(os.path.join(HERE, 'branch.json'), 'w'), indent=1)
