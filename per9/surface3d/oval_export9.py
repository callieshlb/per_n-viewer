#!/usr/bin/env python3
"""E95: real ovals for the viewer (E94 overlay_export.py, ovals only; G_9 edges in overlay_data.json are kept). Only the ovals
listed in OVALS (those that PASSED U6) are exported. For each run of an oval (consecutive lifted points in one tile; oval_lifts.npz
of oval_check9.py) the piece is the holonomy axis (end points xi1, xi2 in that tile's frame) clipped to the Dirichlet polygon D:
in the Klein model the axis is the chord xi1 xi2 and D is an intersection of half-planes (fuchs.halfplanes), so the clip is exact.
Raw checks printed: max distance of the run's points to its axis; sum of the hyperbolic lengths of the distinct pieces against the
geodesic length 2 arccosh(|tr h|/2) of oval_check.json (the pieces must tile the closed geodesic once).
usage: OVALS=4 oval_export9.py [lifts.npz] [oval_check.json]"""
import json, os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from fuchs import to_klein, to_poinc, hdist
from meshloc import load_D
from overlay_export_lib import geod_dist
LIFT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'oval_lifts.npz')
CHK = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, 'oval_check.json')
ids = [int(v) for v in os.environ.get('OVALS', '').split(',') if v != '']
res, D = load_D(here); U, T, _ = fu.halfplanes(res['face_g'])
L = np.load(LIFT); R = {r['oval']: r for r in json.load(open(CHK))}
out = []
for i in ids:
    f, x1, x2, run = L[f'f{i}'], L[f'x1_{i}'], L[f'x2_{i}'], L[f'run{i}']
    pieces = []; worst = 0.0
    for r in np.unique(run):
        idx = np.where(run == r)[0]; a, b = complex(x1[idx[0]]), complex(x2[idx[0]])
        assert np.isfinite(a) and np.isfinite(b)
        t0, t1 = 0.0, 1.0                                    # clip the chord a + t (b - a) to { Re(k conj U_j) <= T_j }
        for u, tt in zip(U, T):
            A = (a*np.conj(u)).real - tt; B = ((b - a)*np.conj(u)).real
            if abs(B) < 1e-300:
                if A > 0: t0, t1 = 1.0, 0.0
                continue
            ts = -A/B
            if B > 0: t1 = min(t1, ts)
            else: t0 = max(t0, ts)
        assert t1 > t0, (i, r, t0, t1)
        P, Q = to_poinc(np.array([a + t0*(b - a), a + t1*(b - a)]))
        dmax = float(geod_dist(f[idx], a, b).max()); worst = max(worst, dmax)
        pieces.append(dict(P=[P.real, P.imag], Q=[Q.real, Q.imag], xi=[[a.real, a.imag], [b.real, b.imag]], n=int(len(idx)), dmax=dmax,
                           length=float(hdist(P, np.array([Q]))[0])))
    uniq = []
    for pc in pieces:
        a = complex(*pc['P']); b = complex(*pc['Q'])
        if any(min(abs(a - complex(*u['P'])) + abs(b - complex(*u['Q'])), abs(a - complex(*u['Q'])) + abs(b - complex(*u['P']))) < 1e-6 for u in uniq): continue
        uniq.append(pc)
    tot = sum(p['length'] for p in uniq); Lg = R[i]['geodesic_length']
    print(f'oval {i}: {len(f)} points, {len(np.unique(run))} runs, {len(uniq)} distinct pieces; max distance of the points to their piece {worst:.3e}; '
          f'piece lengths {[round(p["length"], 6) for p in uniq]}, sum {tot:.9f}, geodesic length {Lg:.9f}, difference {tot - Lg:.3e}')
    out.append(dict(id=i, pieces=uniq, maxdev=worst, geodesic_length=Lg, pieces_length=tot,
                    U6=dict(max_dist_to_axis=R[i]['max_dist_to_axis'], length_ratio_minus_1=R[i]['length_ratio_minus_1'])))
O = json.load(open(os.path.join(here, 'overlay_data.json'))); O['ovals'] = out
json.dump(O, open(os.path.join(here, 'overlay_data.json'), 'w'))
print(f'wrote overlay_data.json: {len(out)} oval(s) {ids}, {len(O["edges"])} edges kept')
