#!/usr/bin/env python3
"""E95 full run: the 252 period-9 centres in the chart (x, y) = (a_8, a_3), from E90's census stage 1 (pern.exact_period_angles,
spider, newton_center at 60 digits; centres deduplicated as census.py) and pern.center_to_chart / pern.orbit as E94's atlas_points.py.
A2 raw checks: |P(x, y)| / sum |terms| (512 bits) and nearest-root distance / runner-up (flint certified roots at y).
Writes explorer/data.json {'centers': [{x, y, ks, real}]} (the input of u5_compare.py). Raw values only."""
import json, os, sys, time
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
E90 = os.path.join(HERE, '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import pern
import trackmp as tm
import options_probe_lib as L
from gmpy2 import mpc
n = 9; dps = 60; mp.mp.dps = dps
t0 = time.time()
centers = []
for k in pern.exact_period_angles(n):
    c, it = pern.spider(k, n)
    cc = pern.newton_center(c, n, dps=dps)
    for e in centers:
        if abs(e['c'] - cc) < mp.mpf(10)**(-dps//2):
            e['ks'].append(k); break
    else:
        centers.append(dict(c=cc, ks=[k]))
print('centres', len(centers), '| angle multiplicities', sorted(set(len(e['ks']) for e in centers)), f'| {time.time()-t0:.1f} s', flush=True)
C = tm.Curve()
out = []; relP = []; ratios = []
mp.mp.dps = 40
for e in centers:
    b, c, _ = pern.center_to_chart(e['c'])
    a = pern.orbit(b, c, n)
    x, y = complex(a[n-1]), complex(a[3])
    cc = C.cy(mpc(y))[0]; X = mpc(x)
    val = sum(ck*X**k for k, ck in enumerate(cc))
    tot = sum(abs(q)*abs(mpc(y))**j*abs(X)**k for k, row in enumerate(C.C) for j, q in enumerate(row))
    relP.append(float(abs(val)/tot))
    rts = sorted(abs(complex(r.mid()) - x) for r in L.froots(C, y))
    ratios.append(rts[1]/max(rts[0], 1e-300))
    out.append(dict(x=[x.real, x.imag], y=[y.real, y.imag], ks=sorted(e['ks']), real=abs(complex(e['c']).imag) < 1e-25))
print(f'A2 raw: max |P|/sum|terms| {max(relP):.2e} | min runner-up/nearest root ratio {min(ratios):.2e} | real centres {sum(o["real"] for o in out)} | min |y| {min(abs(complex(*o["y"])) for o in out):.4f} | {time.time()-t0:.1f} s', flush=True)
os.makedirs(os.path.join(HERE, 'explorer'), exist_ok=True)
json.dump(dict(centers=out), open(os.path.join(HERE, 'explorer', 'data.json'), 'w'))
print('wrote explorer/data.json')
