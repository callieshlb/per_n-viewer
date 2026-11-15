#!/usr/bin/env python3
"""E92-U6: are the 2 real ovals of Per_6(0) geodesics of the computed hyperbolic metric?
Run in a mesh directory after uniformize.py and develop.py (needs surface_mesh.json, sphere_mesh.json, edge_perms.json,
develop.npz). The ovals come from E90 (data_n6_eps1e-4_dps60.json): each is a cyclic chain of traced real arcs between
real punctures; points are re-traced with pern.trace_real (record=...) and mapped to (y, x) = (a_3, a_5) by pern.orbit,
cached in ../oval_points.json (or ../../ for the coarse mesh). Each point is located in the mesh (sphere face by
gnomonic barycentrics, sheet by the closest x), placed in the developed disk (Klein-barycentric in the face), and the
oval is lifted continuously (side pairings applied when it crosses the cut). Outputs per oval: holonomy trace, geodesic
length 2 arccosh(|tr|/2), polyline length, max hyperbolic distance from the lifted oval to the holonomy axis."""
import base64
import json
import os
import sys
import numpy as np
from scipy.spatial import cKDTree
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from develop import Mq, mob_ap, mob_norm, hdist, to_klein, to_poinc
E90 = os.path.join(here, '..', '..', 'inputs', 'E90-real-circles')
CACHE = os.path.join(here, 'oval_points.json')
def trace_ovals():
    sys.path.insert(0, E90)
    import mpmath as mp
    import pern
    d = json.load(open(os.path.join(E90, 'data_n6_eps1e-4_dps60.json')))
    n, dps, eps = d['n'], int(d['dps']), mp.mpf(d['eps'])
    cat, arcs = d['catalog'], d['arcs']
    at = {}
    for ci, C in enumerate(cat): at.setdefault(C['p'], []).append(ci)
    seen, ovals = set(), []
    for c0 in range(len(cat)):
        if c0 in seen: continue
        seq, c = [], c0
        while True:
            e = arcs[c]['end']; seen |= {c, e}; seq.append((c, e))
            c = [k for k in at[cat[e]['p']] if k != e][0]
            if c == c0: break
        ovals.append(seq)
    out = []
    for seq in ovals:
        pts = []
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
            print(f'  arc {c} -> {e}: status {r["status"]}, points so far {len(pts)}')
        out.append(dict(arcs=seq, punctures=[cat[c]['p'] for c, _ in seq], pts=pts))
    json.dump(out, open(CACHE, 'w'))
    return out
if __name__ == '__main__':
    ovals = json.load(open(CACHE)) if os.path.exists(CACHE) else trace_ovals()
    print('ovals:', [(len(o['arcs']), len(o['pts'])) for o in ovals], '(arcs, points)')
    from meshloc import Locator, chord as chord_
    L = Locator('.')
    D = np.load('develop.npz'); z = D['z']; Tn = D['Tn']; G = list(D['gens']); G2 = G + [np.linalg.inv(g) for g in G]
    def candidates(y, x):
        """lifts ranked by exact continuation (root followed from a corner to y inside the triangle); keep ties."""
        r = L.locate(y, x)
        if r is None: return []
        best = r[0]
        return [c for c in L.last_cands if c[0] <= max(10*best, 1e-7)]
    def place_(f, wf):
        a_, b_, c_ = z[Tn[f]]
        Mi = np.linalg.inv(Mq(a_)); Kb, Kc = to_klein(mob_ap(Mi, b_)), to_klein(mob_ap(Mi, c_))
        return mob_ap(Mq(a_), to_poinc(wf[1]*Kb + wf[2]*Kc))
    gens_try = [np.eye(2, dtype=complex)] + G2
    results = []; LIFTS = []
    for oi, ov in enumerate(ovals):
        pts = np.array(ov['pts']); lifts, xmis, h = [], [], np.eye(2, dtype=complex)
        prev, jumps, maxstep = None, 0, 0.0
        # start well away from every puncture (|y|, |x| < 2) so the first sheet choice is unambiguous; continuity resolves
        # the ties met later (near punctures two sheets can arrive at the same x)
        k0 = next(k for k in range(len(pts)) if abs(pts[k, 0]) < 2 and abs(pts[k, 1]) < 2 and len(candidates(complex(pts[k, 0]), complex(pts[k, 1]))) == 1)
        pts = np.concatenate([pts[k0:], pts[:k0]])
        print(f'   oval {oi}: starting at point {k0} (y = {pts[0, 0]:.4f}, x = {pts[0, 1]:.4f})')
        P = list(pts) + [pts[0]]
        amb = 0
        for k, (yr, xr) in enumerate(P):
            cands = candidates(complex(yr), complex(xr))
            if not cands: continue
            if len(cands) > 1: amb += 1
            best = None
            for dx, f, wf in cands:
                w = place_(f, wf)
                if prev is None:
                    cand_lift, cand_h, dist = w, h, 0.0
                else:
                    ds = [float(np.nan_to_num(hdist(prev, mob_ap(h @ s, w)), nan=np.inf)) for s in gens_try]
                    j = int(np.argmin(ds)); cand_h = h @ gens_try[j]; dist = ds[j]
                    if dist > 0.2:
                        d2 = min((hdist(prev, mob_ap(h @ s @ t, w)), si, ti) for si, s in enumerate(G2) for ti, t in enumerate(G2))
                        if d2[0] < dist: cand_h, dist = h @ G2[d2[1]] @ G2[d2[2]], d2[0]
                    cand_lift = mob_ap(cand_h, w)
                if best is None or dist < best[0]: best = (dist, cand_lift, cand_h, dx)
            dist, cur, h, dx = best; xmis.append(dx)
            if prev is not None and (dist > 0.05 or not np.isfinite(dist) or np.abs(h).max() > 1e6):
                print(f'   oval {oi} step {k}: y = {yr:.6g}, x = {xr:.6g}, candidates {len(cands)}, step {dist:.3e}, |h| {np.abs(h).max():.3e}, |lift| {abs(cur):.12f}')
            if prev is not None:
                maxstep = max(maxstep, float(dist))
            lifts.append(cur); prev = cur
        lifts = np.array(lifts)
        g = mob_norm(h); tr = g[0, 0] + g[1, 1]
        Lg = 2*np.arccosh(abs(tr.real)/2) if abs(tr.real) > 2 else 0.0
        Lp = float(np.sum(hdist(lifts[:-1], lifts[1:])))
        # axis: fixed points of g on the unit circle
        a, b, c, dd = g[0, 0], g[0, 1], g[1, 0], g[1, 1]
        xi = np.roots([c, dd - a, -b]); al, be = np.angle(xi[0]), np.angle(xi[1])
        half = (be - al)/2
        if abs(np.cos(half)) < 1e-9:
            th = al; dist = np.arcsinh(2*np.abs((lifts*np.exp(-1j*th)).imag)/(1 - np.abs(lifts)**2))
        else:
            Cc = np.exp(1j*(al + be)/2)/np.cos(half); rho = abs(np.tan(half))
            dist = np.arcsinh(np.abs(np.abs(lifts - Cc)**2 - rho**2)/(rho*(1 - np.abs(lifts)**2)))
        res = dict(oval=oi, arcs=len(ov['arcs']), points=len(lifts), ambiguous_points=amb, max_lift_step=maxstep,
                   max_x_mismatch=float(np.max(xmis)), trace=float(tr.real), trace_imag=float(abs(tr.imag)),
                   geodesic_length=float(Lg), polyline_length=Lp, length_ratio_minus_1=Lp/Lg - 1 if Lg else None,
                   axis_fixed_point_moduli=[float(abs(x)) for x in xi], max_dist_to_axis=float(dist.max()),
                   median_dist_to_axis=float(np.median(dist)))
        results.append(res); LIFTS.append(dict(lifts=lifts, g=g))
        print('U6 oval', oi, json.dumps(res))
    json.dump(results, open('oval_check.json', 'w'), indent=1)
    np.savez('oval_lifts.npz', **{f'lifts{i}': L['lifts'] for i, L in enumerate(LIFTS)}, **{f'g{i}': L['g'] for i, L in enumerate(LIFTS)})
