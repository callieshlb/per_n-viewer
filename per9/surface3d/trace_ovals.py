#!/usr/bin/env python3
"""E94 (n = 8) copy. E93-U6 input (E92 oval_check.trace_ovals with n = 7): the 4 real ovals of Per_7(0) from E90 (data_n8_eps1e-3_dps60.json):
each is a cyclic chain of traced real arcs between real punctures; points re-traced with pern.trace_real (record=...) and mapped
to (y, x) = (a_3, a_6) by pern.orbit. Writes oval_points.json."""
import json
import os
import sys
import time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
E90 = os.path.join(here, '..', '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import mpmath as mp
import pern
d = json.load(open(os.path.join(E90, 'data_n9_eps1e-3_dps60.json')))   # E95: n = 9
n, dps, eps = d['n'], int(d['dps']), mp.mpf(d['eps'])
cat, arcs = d['catalog'], d['arcs']
at = {}
for ci, C in enumerate(cat): at.setdefault(C['p'], []).append(ci)
seen, ovals, dropped = set(), [], []
for c0 in range(len(cat)):
    if c0 in seen: continue
    seq, c = [], c0
    good = True
    while True:
        e = arcs[c]['end']
        if e is None or c in seen and c != c0:          # E95: unlanded arc (maxsteps) or a walk that does not close -> not an oval
            good = False; seen.add(c); break
        seen |= {c, e}; seq.append((c, e))
        oth = [k for k in at[cat[e]['p']] if k != e]
        if not oth: good = False; break
        c = oth[0]
        if c == c0: break
    if good: ovals.append(seq)
    else: dropped.append(len(seq))
print('ovals (arcs per oval):', [len(s) for s in ovals], '| dropped chains (arcs before the break):', dropped, '| E90 circles sizes', [len(c)//2 if isinstance(c, list) else c for c in d['circles']], flush=True)
out = []
t0 = time.time()
for seq in ovals:
    pts = []; arc_end = []
    for c, e in seq:
        rec = []; mp.mp.dps = dps
        r = pern.trace_real(mp.mpf(cat[c]['b']), mp.mpf(cat[c]['c']), n, eps, dps, record=rec)
        for b_, c_ in rec:
            try:
                a = pern.orbit(mp.mpf(b_), mp.mpf(c_), n)
                x, y = complex(a[n - 1]), complex(a[3])
            except ZeroDivisionError:
                continue
            if np.isfinite(x) and np.isfinite(y): pts.append([y.real, x.real])
        arc_end.append(len(pts))
        print(f'  arc {c} -> {e}: status {r["status"]}, points so far {len(pts)} | {time.time()-t0:.0f} s', flush=True)
    out.append(dict(arcs=seq, punctures=[cat[c]['p'] for c, _ in seq], pts=pts, arc_end=arc_end))   # E95: arc_end[i] = number of points up to the end of arc i
json.dump(out, open(os.path.join(here, 'oval_points.json'), 'w'))
print('wrote oval_points.json', [(len(o['arcs']), len(o['pts'])) for o in out])
