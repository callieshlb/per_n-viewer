#!/usr/bin/env python3
"""E92 step B: the 5 charts (sheets of y) of Per_6(0) with the dynamics of the free critical point.
usage: atlas_render.py NPIX [S]
Display coordinate u in the unit disk, y = S*u/(1-|u|) (boundary circle = y = infinity).
Sheet labels: roots x(y) continued from the base point y0 along the straight segment y0 -> y
(cut plane star-shaped about y0; cuts are rays from each branch point away from y0).
Colour: capture of omega by the marked cycle (hue = which marked point it lands next to, a_0..a_5),
type D (grey: period 1 dark = polynomial components; period >= 2 lighter), undecided black.
Writes atlas_{NPIX}.npz and atlas_{NPIX}.png; prints calibration A2/A3 raw values.
"""
import os
import json
import math
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../inputs/E90-real-circles')
import pern  # noqa: E402
import mpmath as mp  # noqa: E402
import sympy as sp  # noqa: E402

NP = int(sys.argv[1]) if len(sys.argv) > 1 else 160
S = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
from atlas_lib import *  # noqa: F401,F403
from atlas_lib import y0, x0, bps, coef, d, roots_batch, match, label

# ---- grid
u1 = np.linspace(-0.999, 0.999, NP)
U = u1[None, :] + 1j*u1[::-1, None]
inside = np.abs(U) < 0.999
Ys = np.where(inside, S*U/(1-np.abs(U)), 0)
flat = Ys[inside]
X, nbad = label(flat)
print(f'pixels {flat.size}, ambiguous matching steps resolved by brute force: {nbad}')

# ---- dynamics
N_IT, DELTA = 600, 1e-4


def chord(a, b):
    return np.abs(a-b)/np.sqrt((1+np.abs(a)**2)*(1+np.abs(b)**2))


def fmap(z, x, w):
    big = np.abs(z) > 1e150
    zs = np.where(big, 1.0, z)
    small = np.abs(zs) < 1e-150
    zs = np.where(small, 1e-150, zs)
    out = 1 - (x+w)/zs + x*w/(zs*zs)
    return np.where(big, 1.0+0j, out)


cls = np.zeros((d,)+flat.shape, dtype=np.int16)   # 1..6 capture at a_{j-1}; 10+p type D period p; 0 undecided
ctime = np.zeros((d,)+flat.shape, dtype=np.int16)
for k in range(d):
    x = X[:, k]
    with np.errstate(all='ignore'):
        w = (x+flat-1)/(x-1)
        om = 2*x*w/(x+w)
        a4 = fmap(flat, x, w)
        marked = [np.zeros_like(x), None, np.ones_like(x), flat, a4, x]
        z = om.copy()
        done = np.zeros(flat.shape, bool)
        for it in range(N_IT):
            z = fmap(z, x, w)
            for j, a in enumerate(marked):
                dj = (1/np.sqrt(1+np.abs(z)**2)) if a is None else chord(z, a)
                hit = (~done) & (dj < DELTA)
                cls[k][hit] = j+1
                ctime[k][hit] = it+1
                done |= hit
        tail = [z]
        for it in range(24):
            z = fmap(z, x, w)
            tail.append(z)
        for p in range(1, 13):
            hit = (~done) & (chord(tail[-1], tail[-1-p]) < 1e-7)
            cls[k][hit] = 10+p
            done |= hit
np.savez_compressed(f'atlas_{NP}.npz', Ys=Ys, inside=inside, X=X, cls=cls, ctime=ctime, S=S)
# A5: labels stable under doubling the continuation steps (sample of pixels)
rng = np.random.default_rng(0)
smp = rng.choice(flat.size, size=min(400, flat.size), replace=False)
X2, _ = label(flat[smp], K1=1400, K2=240)
print('A5 label stability on', len(smp), 'pixels: max |x(K) - x(2K)| per sheet', np.abs(X2 - X[smp]).max(0).round(12).tolist())

# ---- calibration A2: centers
mp.mp.dps = 40
cen = json.load(open(os.path.dirname(os.path.abspath(__file__)) + '/../inputs/E90-real-circles/data_n6_eps1e-4_dps60.json'))
cpts = []
for e in cen['centers']:
    b, c, _ = pern.center_to_chart(mp.mpc(mp.mpf(e['c'][0]), mp.mpf(e['c'][1])))
    a = pern.orbit(b, c, 6)
    cpts.append((complex(a[5]), complex(a[3]), min(e['ks'])))
Xc, _ = label([c[1] for c in cpts])
sheet_c = []
for (xc, yc, k), row in zip(cpts, Xc):
    dist = np.abs(row - xc)
    sheet_c.append(int(dist.argmin()))
    Pval = abs(complex(sum(complex(cf(yc))*xc**(d-i) for i, cf in enumerate(coef))))
print('A2 centers: max |P7(x_c,y_c)| (unnormalised)', max(abs(complex(sum(complex(cf(yc))*xc**(d-i) for i, cf in enumerate(coef)))) for xc, yc, k in cpts),
      '| max over centers of distance to nearest labelled root', max(float(np.abs(r - c[0]).min()) for r, c in zip(Xc, cpts)),
      '| centers per sheet', [sheet_c.count(s) for s in range(d)])
json.dump(dict(centers=[dict(x=[c[0].real, c[0].imag], y=[c[1].real, c[1].imag], k=c[2], sheet=s)
                        for c, s in zip(cpts, sheet_c)]), open('atlas_centers.json', 'w'), indent=1)
hist = {int(v): int((cls == v).sum()) for v in np.unique(cls)}
print('class histogram over all sheets (0 undecided, 1-6 capture at a_0..a_5, 10+p type D period p):', hist)
